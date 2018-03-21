import time
import datetime
import threading
import sys
import subprocess
import logging
import json
import gzip
import os
import traceback

NCBITOOL_S3_PATH = "s3://czbiohub-infectious-disease/ncbitool" # S3 location of ncbitool executable
ACCESSION2TAXID = 's3://czbiohub-infectious-disease/references/accession2taxid.db.gz'
LINEAGE_SHELF = 's3://czbiohub-infectious-disease/references/taxid-lineages.db'
VERSION_NONE = -1

# data directories
ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + '/idseq/data' # generated data go here
REF_DIR = ROOT_DIR + '/idseq/ref' # referene genome / ref databases go here

OUTPUT_VERSIONS = []

print_lock = threading.RLock()

# peak network & storage perf for a typical small instance is saturated by just a few concurrent streams
MAX_CONCURRENT_COPY_OPERATIONS = 8
iostream = threading.Semaphore(MAX_CONCURRENT_COPY_OPERATIONS)


class StatsFile(object):

    def __init__(self, stats_filename, local_results_dir, s3_input_dir, s3_output_dir):
        self.local_results_dir = local_results_dir
        self.stats_filename = stats_filename
        self.s3_input_dir = s3_input_dir
        self.s3_output_dir = s3_output_dir
        self.stats_path = os.path.join(self.local_results_dir, self.stats_filename)
        self.data = []
        self.mutex = threading.RLock()

    def load_from_s3(self):
        stats_s3_path = os.path.join(self.s3_input_dir, self.stats_filename)
        if check_s3_file_presence(stats_s3_path):
            execute_command("aws s3 cp --quiet %s %s/" % (stats_s3_path, self.local_results_dir))
        self._load()

    def _load(self):
        if os.path.isfile(self.stats_path):
            with open(self.stats_path) as f:
                self.data = json.load(f)

    def _save(self):
        with open(self.stats_path, 'wb') as f:
            json.dump(self.data, f)

    def save_to_s3(self):
        with self.mutex:
            self._save()
            execute_command("aws s3 cp --quiet %s %s/" % (self.stats_path, self.s3_output_dir))

    def get_total_reads(self):
        for item in self.data:
            if "total_reads" in item:
                return item["total_reads"]
        # check run star if total reads not available
        for item in self.data:
            if item.get("task") == 'run_star':
                return item["reads_before"]
        # if no entry for run_star, host-filtering was skipped: use "remaining_reads" instead
        return self.get_remaining_reads()
        # previous fall-back value didn't fit into integer type of SQL dfatabase:
        # return 0.1

    def get_remaining_reads(self):
        return (item for item in self.data if item.get("task") == "run_gsnapl_remotely").next().get("reads_before")

    def count_reads(self, func_name, before_filename, before_filetype, after_filename, after_filetype):
        records_before = count_reads(before_filename, before_filetype)
        records_after = count_reads(after_filename, after_filetype)
        new_data = [datum for datum in self.data if datum.get('task') != func_name]
        if len(new_data) != len(self.data):
            write_to_log("Overwriting counts for {}".format(func_name), warning=True)
            self.data = new_data
        self.data.append({'task': func_name, 'reads_before': records_before, 'reads_after': records_after})


class MyThread(threading.Thread):
    def __init__(self, target, args):
        super(MyThread, self).__init__()
        self.args = args
        self.target = target
        self.exception = None
        self.completed = False
        self.result = None

    def run(self):
        try:
            self.result = self.target(*self.args)
            self.exception = False
        except:
            traceback.print_exc()
            self.exception = True
        finally:
            self.completed = True


class Updater(object):

    def __init__(self, update_period, update_function):
        self.update_period = update_period
        self.update_function = update_function
        self.timer_thread = None
        self.t_start = time.time()

    def relaunch(self, initial_launch=False):
        if self.timer_thread and not initial_launch:
            t_elapsed = time.time() - self.t_start
            self.update_function(t_elapsed)
        self.timer_thread = threading.Timer(self.update_period, self.relaunch)
        self.timer_thread.start()

    def __enter__(self):
        self.relaunch(initial_launch=True)
        return self

    def __exit__(self, *args):
        self.timer_thread.cancel()


class CommandTracker(Updater):

    lock = threading.RLock()
    count = 0

    def __init__(self, update_period=15):
        super(CommandTracker, self).__init__(update_period, self.print_update_and_enforce_timeout)
        # User can set the watchdog to a function that takes self.id and t_elapsed as single arg
        self.proc = None
        self.timeout = None
        self.t_sigterm_sent = None
        self.t_sigkill_sent = None
        self.grace_period = update_period / 2.0
        with CommandTracker.lock:
            self.id = CommandTracker.count
            CommandTracker.count += 1

    def print_update_and_enforce_timeout(self, t_elapsed):
        with print_lock:
            if self.proc == None or self.proc.poll() == None:
                print "Command %d still running after %3.1f seconds." % (self.id, t_elapsed)
            else:
                # This should be uncommon, unless there is lengthy python processing following
                # the command in the same CommandTracker "with" block.
                print "Command %d postprocessing after %3.1f seconds." % (self.id, t_elapsed)
            sys.stdout.flush()
        self.enforce_timeout(t_elapsed)

    def enforce_timeout(self, t_elapsed):
        if not self.proc or t_elapsed <= self.timeout or self.proc.poll() != None:
            # Unregistered subprocess, subprocess not yet timed out, or subprocess already exited.
            pass
        elif not self.t_sigterm_sent:
            # Send SIGTERM first.
            with print_lock:
                print "Command %d has exceeded timeout of %3.1f seconds.  Sending SIGTERM." % (self.id, self.timeout)
                sys.stdout.flush()
            self.t_sigterm_sent = time.time()
            self.proc.terminate()
        elif not self.t_sigkill_sent:
            # Grace_period after SIGTERM, send SIGKILL.
            if time.time() > self.t_sigterm_sent + self.grace_period:
                with print_lock:
                    print "Command %d still alive %3.1f seconds after SIGTERM.  Sending SIGKILL." % (self.id, time.time() - self.t_sigterm_sent)
                    sys.stdout.flush()
                self.t_sigkill_sent = time.time()
                self.proc.kill()
        else:
            with print_lock:
                print "Command %d still alive %3.1f seconds after SIGKILL." % (self.id, time.time() - self.t_sigkill_sent)
                sys.stdout.flush()


class ProgressFile(object):

    def __init__(self, progress_file):
        self.progress_file = progress_file
        self.tail_subproc = None

    def __enter__(self):
        # TODO:  Do something else here. Tail gets confused if the file pre-exists.  Also need to rate-limit.
        if self.progress_file:
            self.tail_subproc = subprocess.Popen("touch {pf} ; tail -f {pf}".format(pf=self.progress_file), shell=True)
        return self

    def __exit__(self, *args):
        if self.tail_subproc:
            self.tail_subproc.kill()


def execute_command_with_output(command, progress_file=None, timeout=None, grace_period=None):
    with CommandTracker() as ct:
        with print_lock:
            print "Command {}: {}".format(ct.id, command)
        with ProgressFile(progress_file):
            # Capture only stdout.   Send stderr to parent stderr.   Take input from parent stdin.
            ct.proc = subprocess.Popen(command, shell=True, stdin=sys.stdin.fileno(), stdout=subprocess.PIPE, stderr=sys.stderr.fileno())
            if timeout:
                ct.timeout = timeout
            if grace_period:
                ct.grace_period = grace_period
            stdout, _undefined = ct.proc.communicate()
            if ct.proc.returncode:
                raise subprocess.CalledProcessError(ct.proc.returncode, command, stdout)
            return stdout


def execute_command_realtime_stdout(command, progress_file=None):
    with CommandTracker() as ct:
        with print_lock:
            print "Command {}: {}".format(ct.id, command)
        with ProgressFile(progress_file):
            subprocess.check_call(command, shell=True)


def execute_command(command, progress_file=None):
    execute_command_realtime_stdout(command, progress_file)


def remote_command(base_command, key_path, remote_username, instance_ip):
    return 'ssh -o "StrictHostKeyChecking no" -o "ConnectTimeout 15" -i %s %s@%s "%s"' % (key_path, remote_username, instance_ip, base_command)


def check_s3_file_presence(s3_path):
    command = "aws s3 ls %s | wc -l" % s3_path
    try:
        return int(execute_command_with_output(command).rstrip())
    except:
        return 0


def scp(key_path, remote_username, instance_ip, remote_path, local_path):
    assert " " not in key_path
    assert " " not in remote_path
    assert " " not in local_path
    return 'scp -o "StrictHostKeyChecking no" -o "ConnectTimeout 15" -i {key_path} {username}@{ip}:{remote_path} {local_path}'.format(
        key_path=key_path,
        username=remote_username,
        ip=instance_ip,
        remote_path=remote_path,
        local_path=local_path)


class TimeFilter(logging.Filter):

    def __init__(self, *args, **kwargs):
        super(TimeFilter, self).__init__(*args, **kwargs)
        self.last = None
        self.lock = threading.RLock()

    def filter(self, record):
        with self.lock:
            last = self.last
            if last == None:
                last = record.relativeCreated
            delta = datetime.datetime.fromtimestamp(record.relativeCreated/1000.0) - datetime.datetime.fromtimestamp(last/1000.0)
            record.time_since_last = '{0:.2f}'.format(delta.seconds + delta.microseconds/1000000.0)
            self.last = record.relativeCreated
            return True

def percent_str(percent):
    try:
        return "%3.1f" % percent
    except:
        return str(percent)

def extract_m8_readid(read_id_column):
    # Get raw read ID without our annotations:
    # If m8 is from gsnap (case 1), 1 field has been prepended (taxid field).
    # If m8 is from rapsearch in an old version of the pipeline (case 2), 3 fields have been prepended (taxid, 'NT', NT accession ID).
    # If m8 is from rapsearch in a newer version of the pipeline (case 3), many fields have been prepended (taxid, alignment info fields),
    # but the delimiter ":read_id:" marks the beginning of the raw read ID.
    read_id_column = read_id_column.split(":", 1)[1] # remove taxid field (all cases)
    if ":read_id:" in read_id_column: # case 3
        raw_read_id = read_id_column.split(":read_id:")[1]
    elif read_id_column.startswith("NT:"): # case 2
        raw_read_id = read_id_column.split(":", 2)[2]
    else: # case 1
        raw_read_id = read_id_column
    return raw_read_id

def count_reads(file_name, file_type):
    count = 0
    if file_name[-3:] == '.gz':
        f = gzip.open(file_name)
    else:
        f = open(file_name)
    for line in f:
        if file_type == "fastq_paired":
            count += 2./4
        elif file_type == "fastq":
            count += 1./4
        elif file_type == "fasta_paired":
            if line.startswith('>'):
                count += 2
        elif file_type == "fasta":
            if line.startswith('>'):
                count += 1
        elif file_type == "m8" and line[0] == '#':
            continue
        else:
            count += 1
    f.close()
    return int(count)


def return_merged_dict(dict1, dict2):
    result = dict1.copy()
    result.update(dict2)
    return result


def configure_logger(log_file):
    LOGGER = logging.getLogger()
    LOGGER.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file)
    formatter = logging.Formatter("%(asctime)s (%(time_since_last)ss elapsed): %(message)s")
    handler.addFilter(TimeFilter())
    handler.setFormatter(formatter)
    LOGGER.addHandler(handler)
    # echo to stdout so they get to cloudwatch
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('(%(time_since_last)ss elapsed): %(message)s')
    handler.setFormatter(formatter)
    LOGGER.addHandler(handler)


def fetch_from_s3(source, destination, auto_unzip, allow_s3mi=False, mutex=threading.RLock(), locks={}):  #pylint: disable=dangerous-default-value
    with mutex:
        if os.path.exists(destination):
            if os.path.isdir(destination):
                destination = os.path.join(destination, os.path.basename(source))
        unzip = auto_unzip and destination.endswith(".gz")
        if unzip:
            destination = destination[:-3]
        abspath = os.path.abspath(destination)
        if abspath not in locks:
            locks[abspath] = threading.RLock()
        destination_lock = locks[abspath]
    with destination_lock:
        if os.path.exists(destination):
            # no need to fetch this file from s3, it has been just produced on this instance
            return destination
        try:
            destdir = os.path.dirname(destination)
            if destdir:
                os.makedirs(destdir)
        except OSError, e:
            # It's okay if the parent directory already exists, but all other errors are fatal.
            if e.errno != os.errno.EEXIST:
                raise
        with iostream:
            try:
                if allow_s3mi:
                    try:
                        install_s3mi()
                    except:
                        allow_s3mi = False
                if unzip:
                    try:
                        assert allow_s3mi
                        execute_command("s3mi cat %s | gzip -dc > %s" % (source, destination))
                    except:
                        execute_command("aws s3 cp --quiet %s - | gzip -dc > %s" % (source, destination))
                else:
                    try:
                        assert allow_s3mi
                        execute_command("s3mi cp %s %s" % (source, destination))
                    except:
                        execute_command("aws s3 cp --quiet %s %s" % (source, destination))
                return destination
            except subprocess.CalledProcessError:
                # Most likely the file doesn't exist in S3.
                return None


def install_s3mi(installed={}, mutex=threading.RLock()): #pylint: disable=dangerous-default-value
    with mutex:
        if installed:
            return
        try:
            # This is typically a no-op.
            execute_command("which s3mi || pip install git+git://github.com/chanzuckerberg/s3mi.git")
            execute_command("s3mi tweak-vm || echo s3mi tweak-vm is impossible under docker")
        finally:
            installed['time'] = time.time()


def fetch_lazy_result(source, destination, allow_s3mi=False):
    return fetch_from_s3(source, destination, auto_unzip=False, allow_s3mi=allow_s3mi) != None


def fetch_reference(source, auto_unzip=True):
    path = fetch_from_s3(source, REF_DIR, auto_unzip, allow_s3mi=True)
    assert path != None
    return path


run_and_log_mutex = threading.RLock()

def run_and_log(logparams, target_outputs, lazy_run, func_name, *args):
    with run_and_log_mutex:
        return run_and_log_work(logparams, target_outputs, lazy_run, func_name, os.path.isfile, *args)

def run_and_log_eager(logparams, target_outputs, func_name, *args):
    with run_and_log_mutex:
        return run_and_log_work(logparams, target_outputs, False, func_name, None, *args)

def run_and_log_s3(logparams, target_outputs, lazy_run, s3_result_dir, func_name, *args):
    def lazy_fetch(output_file):
        return fetch_lazy_result(os.path.join(s3_result_dir, os.path.basename(output_file)), output_file, allow_s3mi=True)
    with run_and_log_mutex:
        return run_and_log_work(logparams, target_outputs, lazy_run, func_name, lazy_fetch, *args)

def run_and_log_work(logparams, target_outputs, lazy_run, func_name, lazy_fetch, *args):
    global OUTPUT_VERSIONS
    global run_and_log_mutex

    LOGGER = logging.getLogger()
    LOGGER.info("========== %s ==========", logparams.get("title"))

    # record version of any reference index used
    version_file_s3 = logparams.get("version_file_s3")
    if version_file_s3:
        version_json = execute_command_with_output("aws s3 cp --quiet %s -" % version_file_s3)
        OUTPUT_VERSIONS.append(json.loads(version_json))

    # upload version output
    output_version_file = logparams.get("output_version_file")
    if output_version_file:
        with open(output_version_file, 'wb') as f:
            json.dump(OUTPUT_VERSIONS, f)
        execute_command("aws s3 cp --quiet %s %s/;" % (output_version_file, logparams["sample_s3_output_path"]))

    # produce the output
    # This is the slow part that happens outside the mutex
    run_and_log_mutex.release()
    try:
        was_lazy = lazy_run and all(lazy_fetch(output) for output in target_outputs)
        if not was_lazy:
            func_name(*args)
    finally:
        run_and_log_mutex.acquire()

    if was_lazy:
        LOGGER.info("output exists, lazy run")
    else:
        LOGGER.info("uploaded output")


def upload_log_file(sample_s3_output_path, lock=threading.RLock()):
    with lock:
        logh = logging.getLogger().handlers[0]
        logh.flush()
        execute_command("aws s3 cp --quiet %s %s/;" % (logh.baseFilename, sample_s3_output_path))


def write_to_log(message, warning=False):
    LOGGER = logging.getLogger()
    with print_lock:
        if warning:
            LOGGER.warn(message)
        else:
            LOGGER.info(message)
        sys.stdout.flush()


def unbuffer_stdout():
    # Unbuffer stdout and redirect stderr into stdout. This helps observe logged events in realtime.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    os.dup2(sys.stdout.fileno(), sys.stderr.fileno())

def upload_commit_sha(version):
    sha_file = os.environ.get('COMMIT_SHA_FILE')
    s3_destination = os.environ.get('OUTPUT_BUCKET')
    if sha_file is None or s3_destination is None:
        return
    sha_file_parts = os.path.splitext(os.path.basename(sha_file))
    aws_batch_job_id = os.environ.get('AWS_BATCH_JOB_ID', 'local')
    sha_file_new_name = "%s_job-%s%s" % (sha_file_parts[0], aws_batch_job_id, sha_file_parts[1])
    execute_command("aws s3 cp --quiet %s %s/%s;" % (sha_file, s3_destination.rstrip('/'), sha_file_new_name))

    # also initialize the main version output file with the commit sha and the job_id
    global OUTPUT_VERSIONS
    aws_batch_job_id = os.environ.get('AWS_BATCH_JOB_ID', 'local')
    OUTPUT_VERSIONS.append({"name": "job_id", "version": aws_batch_job_id})
    with open(sha_file, 'r') as f:
        commit_sha = f.read().rstrip()
    OUTPUT_VERSIONS.append({"name": "idseq-pipeline", "version": version, "commit-sha": commit_sha})

def install_ncbitool_locally(local_work_dir):
    execute_command("aws s3 cp --quiet %s %s/" % (NCBITOOL_S3_PATH, local_work_dir))
    execute_command("chmod u+x %s/ncbitool" % local_work_dir)
    return "%s/ncbitool" % local_work_dir

def install_ncbitool(local_work_dir, remote_work_dir=None, key_path=None, remote_username=None, server_ip=None, sudo=False):
    local_result = install_ncbitool_locally(local_work_dir)
    if remote_work_dir is None:
        return local_result
    sudo_prefix = "sudo " if sudo else ""
    command = sudo_prefix + "aws s3 cp --quiet %s %s/; " % (NCBITOOL_S3_PATH, remote_work_dir)
    command += sudo_prefix + "chmod u+x %s/ncbitool" % (remote_work_dir)
    execute_command(remote_command(command, key_path, remote_username, server_ip))
    remote_result = "%s/ncbitool" % remote_work_dir
    return local_result, remote_result

def get_reference_version_number(ncbitool_path, input_fasta_ncbi_path):
    command = "%s file history %s" % (ncbitool_path, input_fasta_ncbi_path)
    output = execute_command_with_output(command).split("File History: ")[1]
    version_history = json.loads(output)
    version_numbers = [entry["Version"] for entry in version_history]
    return max(version_numbers) # use latest available version

def download_reference_locally(ncbitool_path, input_fasta_ncbi_path, version_number, destination_dir):
    command = "cd %s; %s file --download --version-num %s %s" % (destination_dir, ncbitool_path, version_number, input_fasta_ncbi_path)
    execute_command(command)
    return os.path.join(destination_dir, os.path.basename(input_fasta_ncbi_path))

def download_reference_locally_with_version_any_source_type(ref_file, dest_dir, ncbitool_dest_dir):
    # Get input reference and version number.
    # If download does not use ncbitool (e.g. direct s3 link), indicate that there is no versioning.
    input_fasta_name = os.path.basename(ref_file)
    if ref_file.startswith("s3://"):
        execute_command("aws s3 cp --quiet %s %s/" % (ref_file, dest_dir))
        input_fasta_local = os.path.join(dest_dir, input_fasta_name)
        version_number = VERSION_NONE
    elif ref_file.startswith("ftp://"):
        execute_command("cd %s; wget %s" % (dest_dir, ref_file))
        input_fasta_local = os.path.join(dest_dir, input_fasta_name)
        version_number = VERSION_NONE
    else:
        ncbitool_path = install_ncbitool(ncbitool_dest_dir)
        version_number = get_reference_version_number(ncbitool_path, ref_file)
        input_fasta_local = download_reference_locally(ncbitool_path, ref_file, version_number, dest_dir)
    return input_fasta_local, version_number

def download_reference_on_remote(ncbitool_path, input_fasta_ncbi_path, version_number, destination_dir,
                                 key_path, remote_username, server_ip):
    command = "cd %s; sudo %s file --download --version-num %s %s" % (destination_dir, ncbitool_path, version_number, input_fasta_ncbi_path)
    execute_command(remote_command(command, key_path, remote_username, server_ip))
    return os.path.join(destination_dir, os.path.basename(input_fasta_ncbi_path))

def download_reference_on_remote_with_version_any_source_type(
        ref_file, dest_dir,
        local_ncbitool_dest_dir, remote_ncbitool_dest_dir,
        key_path, remote_username, server_ip, sudo=False):
    # Get input reference and version number.
    # If download does not use ncbitool (e.g. direct s3 link), indicate that there is no versioning.
    input_fasta_name = os.path.basename(ref_file)
    if ref_file.startswith("s3://"):
        input_fasta_remote = os.path.join(dest_dir, input_fasta_name)
        execute_command(remote_command("aws s3 cp --quiet %s %s/" % (ref_file, dest_dir),
                                       key_path, remote_username, server_ip))
        version_number = VERSION_NONE
    elif ref_file.startswith("ftp://"):
        input_fasta_remote = os.path.join(dest_dir, input_fasta_name)
        execute_command(remote_command("cd %s; sudo wget %s" % (dest_dir, ref_file),
                                       key_path, remote_username, server_ip))
        version_number = VERSION_NONE
    else:
        local_ncbitool, remote_ncbitool = install_ncbitool(
            local_ncbitool_dest_dir, remote_ncbitool_dest_dir,
            key_path, remote_username, server_ip, sudo)
        version_number = get_reference_version_number(local_ncbitool, ref_file)
        input_fasta_remote = download_reference_on_remote(
            remote_ncbitool, ref_file, version_number, dest_dir,
            key_path, remote_username, server_ip)
    return input_fasta_remote, version_number

def upload_version_tracker(source_file, output_name, reference_version_number, output_path_s3, indexing_version):
    version_tracker_file = "%s.version.txt" % output_name
    version_json = {"name": output_name,
                    "source_file": source_file,
                    "source_version": reference_version_number,
                    "indexing_version": indexing_version,
                    "generation_date": datetime.datetime.now().isoformat()}
    with open(version_tracker_file, 'wb') as f:
        json.dump(version_json, f)
    execute_command("aws s3 cp --quiet %s %s/" % (version_tracker_file, output_path_s3))
