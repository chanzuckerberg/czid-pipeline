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
import re
import multiprocessing
from functools import wraps
import random

from idseq_pipeline import __version__ as PIPELINE_VERSION
bucket = "s3://idseq-database"
NCBITOOL_S3_PATH = bucket + "/ncbitool"  # S3 location of ncbitool executable
base_dt = '2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800' \
          '-unixtime'
base_s3 = bucket + "/alignment_data"
ACCESSION2TAXID = ("%s/%s/accession2taxid.db" % (base_s3, base_dt))
base_s3 = bucket + "/taxonomy"
LINEAGE_SHELF = ("%s/%s/taxid-lineages.db" % (base_s3, base_dt))

# Don't run into the -2e9 limit of mysql database. Current largest taxid is
# around 2e6 so should be fine.
INVALID_CALL_BASE_ID = -(10**8)

VERSION_NONE = -1

# Data directories
ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + '/idseq/data'  # Generated data go here
REF_DIR = ROOT_DIR + '/idseq/ref'  # Reference genome / ref databases go here

OUTPUT_VERSIONS = []

print_lock = multiprocessing.RLock()

# Peak network and storage perf for a typical small instance is saturated by
# just a few concurrent streams.
MAX_CONCURRENT_COPY_OPERATIONS = 8
iostream = multiprocessing.Semaphore(MAX_CONCURRENT_COPY_OPERATIONS)
# Make a second semaphore for uploads to reserve some capacity for downloads.
MAX_CONCURRENT_UPLOAD_OPERATIONS = 4
iostream_uploads = multiprocessing.Semaphore(MAX_CONCURRENT_UPLOAD_OPERATIONS)

# Definitions for integration with web app
TAX_LEVEL_SPECIES = 1
TAX_LEVEL_GENUS = 2
TAX_LEVEL_FAMILY = 3
NULL_SPECIES_ID = -100
NULL_GENUS_ID = -200
NULL_FAMILY_ID = -300
NULL_LINEAGE = (str(NULL_SPECIES_ID), str(NULL_GENUS_ID), str(NULL_FAMILY_ID))
JOB_SUCCEEDED = "succeeded"
JOB_FAILED = "failed"

AWS_BATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID', 'local')


class StatsFile(object):
    """StatsFile is for handling statistics (e.g. task, reads before,
    reads after) in external command execution.
    """
    def __init__(self, stats_filename, local_results_dir, s3_input_dir,
                 s3_output_dir):
        self.local_results_dir = local_results_dir
        self.stats_filename = stats_filename
        self.s3_input_dir = s3_input_dir
        self.s3_output_dir = s3_output_dir
        self.stats_path = os.path.join(self.local_results_dir,
                                       self.stats_filename)
        self.data = []
        self.mutex = threading.RLock()

    def load_from_s3(self):
        stats_s3_path = os.path.join(self.s3_input_dir, self.stats_filename)
        try:
            cmd = "aws s3 cp --quiet %s %s/" % (stats_s3_path,
                                                self.local_results_dir)
            execute_command(cmd)
        except:
            time.sleep(1.0)
            cmd = "aws s3 cp --quiet %s %s/" % (stats_s3_path,
                                                self.local_results_dir)
            execute_command(cmd)
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
            cmd = "aws s3 cp --quiet %s %s/" % (self.stats_path,
                                                self.s3_output_dir)
            execute_command(cmd)

    def get_item_value(self, key):
        for item in self.data:
            if item.get(key):
                return item[key]

    def get_item_for_task(self, function_name):
        for item in self.data:
            if item.get("task") == function_name:
                return item

    def get_total_reads(self):
        # New style
        tr = self.get_item_value("total_reads")
        if tr is not None:
            return int(tr)
        # For compatibility
        item = self.get_item_for_task("run_star")
        if item is not None:
            return int(item["reads_before"])

    def get_remaining_reads(self):
        return int(self.get_item_value('remaining_reads'))

    def gsnap_ran_in_host_filtering(self):
        return self.get_item_for_task("run_gsnap_filter") is not None

    def count_reads(self, func_name, before_filename, before_filetype,
                    after_filename, after_filetype):
        records_before = count_reads(before_filename, before_filetype)
        records_after = count_reads(after_filename, after_filetype)
        new_data = [
            datum for datum in self.data if datum.get('task') != func_name
        ]

        if len(new_data) != len(self.data):
            msg = "Overwriting counts for {}".format(func_name)
            write_to_log(msg, warning=True)
            self.data = new_data
        self.data.append({
            'task': func_name,
            'reads_before': records_before,
            'reads_after': records_after
        })

    def set_remaining_reads(self):
        last_entry = self.data[-1] or {}
        # We use [] and not get() because, if the last entry doesn't have
        # reads_after, then we have a bug and should crash rather than produce
        # misleading results.
        remaining = last_entry['reads_after']
        self.data.append({'remaining_reads': remaining})


class MyThread(threading.Thread):
    def __init__(self, target, args, kwargs=None):
        super(MyThread, self).__init__()
        self.args = args
        if kwargs is None:
            self.kwargs = {}
        else:
            self.kwargs = kwargs
        self.target = target
        self.exception = None
        self.completed = False
        self.result = None

    def run(self):
        try:
            self.result = self.target(*self.args, **self.kwargs)
            self.exception = False
        except:
            traceback.print_exc()
            self.exception = True
        finally:
            self.completed = True


class Updater(object):
    """Base for CommandTracker."""

    def __init__(self, update_period, update_function):
        self.update_period = update_period
        self.update_function = update_function
        self.timer_thread = None
        self.exited = False
        self.t_start = time.time()

    def relaunch(self, initial_launch=False):
        if self.exited:
            return
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
        self.exited = True


class CommandTracker(Updater):
    """CommandTracker is for running external and remote commands and
    monitoring their progress with log updates and timeouts.
    """
    lock = multiprocessing.RLock()
    count = multiprocessing.Value('i', 0)

    def __init__(self, update_period=15):
        super(CommandTracker, self).__init__(
            update_period, self.print_update_and_enforce_timeout)
        # User can set the watchdog to a function that takes self.id and
        # t_elapsed as single arg
        self.proc = None  # Value indicates registered subprocess.
        self.timeout = None
        self.t_sigterm_sent = None  # First sigterm, then sigkill.
        self.t_sigkill_sent = None
        self.grace_period = update_period / 2.0
        with CommandTracker.lock:
            self.id = CommandTracker.count.value
            CommandTracker.count.value += 1

    def print_update_and_enforce_timeout(self, t_elapsed):
        """Log an update after every polling period to indicate the command is
        still active.
        """
        with print_lock:
            if self.proc is None or self.proc.poll() is None:
                print("Command %d still running after %3.1f seconds." %
                      (self.id, t_elapsed))
            else:
                # This should be uncommon, unless there is lengthy python
                # processing following the command in the same CommandTracker
                # "with" block. Note: Not to be confused with post-processing
                # on the data.
                print("Command %d still postprocessing after %3.1f seconds." %
                      (self.id, t_elapsed))
            sys.stdout.flush()
        self.enforce_timeout(t_elapsed)

    def enforce_timeout(self, t_elapsed):
        """Check the timeout and send SIGTERM then SIGKILL to end a command's
        execution.
        """
        if self.timeout is None or not self.proc or \
                t_elapsed <= self.timeout or self.proc.poll() is not None:
            # Skip if unregistered subprocess, subprocess not yet timed out,
            # or subprocess already exited.
            pass
        elif not self.t_sigterm_sent:
            # Send SIGTERM first.
            with print_lock:
                msg = "Command %d has exceeded timeout of %3.1f seconds. " \
                      "Sending SIGTERM." % (self.id, self.timeout)
                print(msg)
                sys.stdout.flush()
            self.t_sigterm_sent = time.time()
            self.proc.terminate()
        elif not self.t_sigkill_sent:
            # Grace_period after SIGTERM, send SIGKILL.
            if time.time() > self.t_sigterm_sent + self.grace_period:
                with print_lock:
                    msg = "Command %d still alive %3.1f seconds after " \
                          "SIGTERM. Sending SIGKILL." % (self.id, time.time() - self.t_sigterm_sent)
                    print(msg)
                    sys.stdout.flush()
                self.t_sigkill_sent = time.time()
                self.proc.kill()
        else:
            with print_lock:
                msg = "Command %d still alive %3.1f seconds after " \
                      "SIGKILL." % (self.id, time.time() - self.t_sigkill_sent)
                print(msg)
                sys.stdout.flush()


class ProgressFile(object):
    def __init__(self, progress_file):
        self.progress_file = progress_file
        self.tail_subproc = None

    def __enter__(self):
        # TODO: Do something else here. Tail gets confused if the file
        # pre-exists. Also need to rate-limit.
        if self.progress_file:
            self.tail_subproc = subprocess.Popen(
                "touch {pf} ; tail -f {pf}".format(pf=self.progress_file),
                shell=True)
        return self

    def __exit__(self, *args):
        if self.tail_subproc:
            # TODO: Do we need to join the tail subproc after killing it?
            self.tail_subproc.kill()


def run_in_subprocess(target):
    """
    Decorator that executes a function synchronously in a subprocess.
    Use case:

        thread 1:
             compute_something(x1, y1, z1)

        thread 2:
             compute_something(x2, y2, z2)

        thread 3:
             compute_something(x3, y3, z3)

    If compute_something() is CPU-intensive, the above threads won't really run
    in parallel because of the Python global interpreter lock (GIL). To avoid
    this problem without changing the above code, simply decorate the definition
    of compute_something() like so:

         @run_in_subprocess
         def compute_something(x, y, z):
             ...

    Typical subprocess limitations or caveats apply:
       a. The caller can't see any value returned by the decorated function.
          It should output to a file, a pipe, or a multiprocessing queue.
       b. Changes made to global variables won't be seen by parent process.
       c. Use multiprocessing semaphores/locks/etc, not their threading versions.

    Tip: If input from the same file is needed in all invocations of the
    decorated function, do the I/O before the first call, to avoid accessing
    the file multiple times.
    """

    @wraps(target)
    def wrapper(*args, **kwargs):
        p = multiprocessing.Process(target=target, args=args, kwargs=kwargs)
        p.start()
        p.join()
        if p.exitcode != 0:
            raise RuntimeError("Failed {} on {}, {}".format(
                target.__name__, args, kwargs))

    return wrapper


def retry(operation, randgen=random.Random().random):
    """Retry decorator for external commands."""
    # Note the use of a separate random generator for retries so transient
    # errors won't perturb other random streams used in the application.
    @wraps(operation)
    def wrapped_operation(*args, **kwargs):
        remaining_attempts = 3
        delay = 1.0
        while remaining_attempts > 1:
            try:
                return operation(*args, **kwargs)
            except:
                # Random jitter and exponential delay
                time.sleep(delay * (1.0 + randgen()))
                delay *= 3.0
                remaining_attempts -= 1
        # The last attempt is outside try/catch so caller can handle exception
        return operation(*args, **kwargs)

    return wrapped_operation


def execute_command(command,
                    progress_file=None,
                    timeout=None,
                    grace_period=None,
                    capture_stdout=False,
                    merge_stderr=False):
    """Primary way to start external commands in subprocesses and handle
    execution with logging.
    """
    with CommandTracker() as ct:
        with print_lock:
            print("Command {}: {}".format(ct.id, command))
        with ProgressFile(progress_file):
            if timeout:
                ct.timeout = timeout
            if grace_period:
                ct.grace_period = grace_period
            if capture_stdout:
                # Capture only stdout. Child stderr = parent stderr unless
                # merge_stderr specified. Child input = parent stdin.
                ct.proc = subprocess.Popen(
                    command,
                    shell=True,
                    stdin=sys.stdin.fileno(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT
                    if merge_stderr else sys.stderr.fileno())
                stdout, _ = ct.proc.communicate()
            else:
                # Capture nothing. Child inherits parent stdin/out/err.
                ct.proc = subprocess.Popen(command, shell=True)
                ct.proc.wait()
                stdout = None

            if ct.proc.returncode:
                raise subprocess.CalledProcessError(ct.proc.returncode,
                                                    command, stdout)
            if capture_stdout:
                return stdout


def execute_command_with_output(command,
                                progress_file=None,
                                timeout=None,
                                grace_period=None,
                                merge_stderr=False):
    return execute_command(
        command,
        progress_file,
        timeout,
        grace_period,
        capture_stdout=True,
        merge_stderr=merge_stderr)


def remote_command(base_command, key_path, remote_username, instance_ip):
    # ServerAliveInterval to fix issue with containers keeping open an SSH
    # connection even after worker machines had finished running.
    return 'ssh -o "StrictHostKeyChecking no" -o "ConnectTimeout 15" ' \
           '-o "ServerAliveInterval 60" -i %s %s@%s "%s"' % (
        key_path, remote_username, instance_ip, base_command)


def scp(key_path, remote_username, instance_ip, remote_path, local_path):
    assert " " not in key_path
    assert " " not in remote_path
    assert " " not in local_path
    # ServerAliveInterval to fix issue with containers keeping open an SSH
    # connection even after worker machines had finished running.
    return 'scp -o "StrictHostKeyChecking no" -o "ConnectTimeout 15" ' \
           '-o "ServerAliveInterval 60" -i {key_path} ' \
           '{username}@{ip}:{remote_path} {local_path}'.format(
        key_path=key_path,
        username=remote_username,
        ip=instance_ip,
        remote_path=remote_path,
        local_path=local_path)


def percent_str(percent):
    try:
        return "%3.1f" % percent
    except:
        return str(percent)


def count_reads(file_name, file_type, max_reads=None):
    """Count number of reads/records in different file types."""
    count = 0
    if file_name[-3:] == '.gz':
        f = gzip.open(file_name)
    else:
        f = open(file_name)

    for line in f:
        if file_type == "fastq_paired":
            count += 2. / 4
        elif file_type == "fastq":
            count += 1. / 4
        elif file_type == "fasta_paired" and line.startswith('>'):
            count += 2
        elif file_type == "m8" and line[0] == '#':  # Comment lines
            pass
        # elif file_type == "fasta" and line.startswith('>'):
        #     count += 1
        else:
            count += 1

        if max_reads is not None and count > max_reads:
            # If more than max, return early. Using > because of floating point.
            break
    f.close()
    return int(count)


def return_merged_dict(dict1, dict2):
    result = dict1.copy()
    result.update(dict2)
    return result


def configure_logger(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    handler = logging.FileHandler(log_file)
    formatter = logging.Formatter("%(asctime)s: %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # Echo to stdout so they get to CloudWatch
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def fetch_from_s3(src,
                  dst,
                  auto_unzip,
                  allow_s3mi=False,
                  mutex=threading.RLock(),
                  locks={}):  #pylint: disable=dangerous-default-value
    """Fetch a file from S3 if needed using either s3mi or aws cp."""
    with mutex:
        if os.path.exists(dst) and os.path.isdir(dst):
            dst = os.path.join(dst, os.path.basename(src))
        unzip = auto_unzip and dst.endswith(".gz")
        if unzip:
            dst = dst[:-3]  # Remove .gz
        abspath = os.path.abspath(dst)
        if abspath not in locks:
            locks[abspath] = threading.RLock()
        destination_lock = locks[abspath]

    with destination_lock:
        if os.path.exists(dst):
            # No need to fetch this file from s3, it has been just produced
            # on this instance.
            return dst

        try:
            destdir = os.path.dirname(dst)
            if destdir:
                os.makedirs(destdir)
        except OSError as e:
            # It's okay if the parent directory already exists, but all other
            # errors are fatal.
            if e.errno != os.errno.EEXIST:
                print("Error in creating destination directory.")
                raise

        with iostream:
            try:
                if allow_s3mi:
                    try:
                        install_s3mi()
                    except:
                        print("s3mi failed to install.")
                        allow_s3mi = False

                pipe_filter = ""
                if unzip:
                    pipe_filter = "| gzip -dc "
                try:
                    assert allow_s3mi
                    cmd = "s3mi cat {source} {pipe_filter} > {destination}".format(
                        source=src, pipe_filter=pipe_filter, destination=dst)
                    execute_command(cmd)
                except:
                    print(
                        "Failed to download with s3mi. Trying with aws s3 cp..."
                    )
                    cmd = "aws s3 cp --quiet {source} - {pipe_filter} > {destination}".format(
                        source=src, pipe_filter=pipe_filter, destination=dst)
                    execute_command(cmd)
                return dst
            except subprocess.CalledProcessError:
                # Most likely the file doesn't exist in S3.
                print(
                    "Failed to fetch file from S3. Most likely does not exist."
                )
                return None


def install_s3mi(installed={}, mutex=threading.RLock()):  #pylint: disable=dangerous-default-value
    with mutex:
        if installed:  # Mutable default value persists
            return
        try:
            # This is typically a no-op.
            execute_command(
                "which s3mi || pip install git+git://github.com/chanzuckerberg/s3mi.git"
            )
            execute_command(
                "s3mi tweak-vm || echo s3mi tweak-vm is impossible under docker. Continuing..."
            )
        finally:
            installed['time'] = time.time()


def fetch_lazy_result(source, destination, allow_s3mi=False):
    res = fetch_from_s3(
        source, destination, auto_unzip=False,
        allow_s3mi=allow_s3mi) is not None
    return res


def fetch_reference(source, auto_unzip=True, allow_s3mi=True):
    path = fetch_from_s3(source, REF_DIR, auto_unzip, allow_s3mi=allow_s3mi)
    assert path is not None, "Failed to fetch reference from S3."
    return path


run_and_log_mutex = threading.RLock()


def run_and_log(log_params, target_outputs, lazy_run, func_name, *args):
    with run_and_log_mutex:
        return run_and_log_work(log_params, target_outputs, lazy_run,
                                func_name, os.path.isfile, *args)


def run_and_log_eager(log_params, target_outputs, func_name, *args):
    with run_and_log_mutex:
        return run_and_log_work(log_params, target_outputs, False, func_name,
                                None, *args)


def run_and_log_s3(log_params, target_outputs, lazy_run, s3_result_dir,
                   func_name, *args):
    def lazy_fetch(output_file):
        return fetch_lazy_result(
            os.path.join(s3_result_dir, os.path.basename(output_file)),
            output_file,
            allow_s3mi=True)

    with run_and_log_mutex:
        return run_and_log_work(log_params, target_outputs, lazy_run,
                                func_name, lazy_fetch, *args)


def run_and_log_work(log_params, target_outputs, lazy_run, func_name,
                     lazy_fetch, *args):
    global OUTPUT_VERSIONS
    global run_and_log_mutex

    msg = "========== %s ==========" % log_params.get("title")
    write_to_log(msg)

    # record version of any reference index used
    version_file_s3 = log_params.get("version_file_s3")
    if version_file_s3:
        version_json = execute_command_with_output(
            "aws s3 cp --quiet %s -" % version_file_s3)
        OUTPUT_VERSIONS.append(json.loads(version_json))

    # upload version output
    output_version_file = log_params.get("output_version_file")
    if output_version_file:
        with open(output_version_file, 'wb') as f:
            json.dump(OUTPUT_VERSIONS, f)
        execute_command("aws s3 cp --quiet %s %s/;" %
                        (output_version_file,
                         log_params["sample_s3_output_path"]))

    # produce the output
    # This is the slow part that happens outside the mutex
    run_and_log_mutex.release()
    try:
        was_lazy = lazy_run and all(
            lazy_fetch(output) for output in target_outputs)
        if not was_lazy:
            func_name(*args)
    finally:
        run_and_log_mutex.acquire()

    if was_lazy:
        write_to_log("output exists, lazy run")
    else:
        write_to_log(
            "non-lazy run completed, output may be in the process of being uploaded"
        )

    return was_lazy


def upload_log_file(sample_s3_output_path, lock=threading.RLock()):
    with lock:
        logh = logging.getLogger().handlers[0]
        logh.flush()
        execute_command("aws s3 cp --quiet %s %s/;" % (logh.baseFilename,
                                                       sample_s3_output_path))


def write_to_log(message, warning=False, flush=True):
    logger = logging.getLogger()
    with print_lock:
        if warning:
            logger.warning(message)
        else:
            logger.info(message)
        if flush:
            sys.stdout.flush()


def set_up_stdout():
    # Unbuffer stdout and redirect stderr into stdout. This helps observe logged events in realtime.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    os.dup2(sys.stdout.fileno(), sys.stderr.fileno())


def set_up_commit_sha(version, s3_destination=None):
    """Upload a file with the commit hash and also initialize the main version output file.
    """
    sha_file = os.environ.get('COMMIT_SHA_FILE')
    s3_destination = s3_destination or os.environ.get('OUTPUT_BUCKET')
    if sha_file is None or s3_destination is None:
        return
    sha_file_parts = os.path.splitext(os.path.basename(sha_file))
    aws_batch_job_id = AWS_BATCH_JOB_ID
    sha_file_new_name = "%s_job-%s%s" % (sha_file_parts[0], aws_batch_job_id,
                                         sha_file_parts[1])
    execute_command("aws s3 cp --quiet %s %s/%s;" %
                    (sha_file, s3_destination.rstrip('/'), sha_file_new_name))

    # also initialize the main version output file with the commit sha and the job_id
    global OUTPUT_VERSIONS
    OUTPUT_VERSIONS.append({"name": "job_id", "version": aws_batch_job_id})

    with open(sha_file, 'r') as f:
        commit_sha = f.read().rstrip()

    OUTPUT_VERSIONS.append({
        "name": "idseq-pipeline",
        "version": version,
        "commit-sha": commit_sha
    })


def install_ncbitool_locally(local_work_dir):
    if not os.path.isfile(local_work_dir + "/ncbitool"):
        execute_command("aws s3 cp --quiet %s %s/" % (NCBITOOL_S3_PATH,
                                                      local_work_dir))
        execute_command("chmod u+x %s/ncbitool" % local_work_dir)
    return "%s/ncbitool" % local_work_dir


def install_ncbitool(local_work_dir,
                     remote_work_dir=None,
                     key_path=None,
                     remote_username=None,
                     server_ip=None,
                     sudo=False):
    local_result = install_ncbitool_locally(local_work_dir)
    if remote_work_dir is None:
        return local_result
    sudo_prefix = "sudo " if sudo else ""
    command = sudo_prefix + "aws s3 cp --quiet %s %s/; " % (NCBITOOL_S3_PATH,
                                                            remote_work_dir)
    command += sudo_prefix + "chmod u+x %s/ncbitool" % (remote_work_dir)
    execute_command(
        remote_command(command, key_path, remote_username, server_ip))
    remote_result = "%s/ncbitool" % remote_work_dir
    return local_result, remote_result


def get_reference_version_number(ncbitool_path, input_fasta_ncbi_path):
    command = "%s file history %s" % (ncbitool_path, input_fasta_ncbi_path)
    output = execute_command_with_output(command).split("File History: ")[1]
    version_history = json.loads(output)
    version_numbers = [entry["Version"] for entry in version_history]
    return max(version_numbers)  # use latest available version


def download_ref_local(ncbitool_path, input_fasta_ncbi_path,
                       version_number, destination_dir):
    command = "cd %s; %s file --download --version-num %s %s" % (
        destination_dir, ncbitool_path, version_number, input_fasta_ncbi_path)
    execute_command(command)
    return os.path.join(destination_dir,
                        os.path.basename(input_fasta_ncbi_path))


def download_ref_local_with_version_any_type(
        ref_file, dest_dir, ncbitool_dest_dir, auto_unzip=False):
    # Get input reference and version number. If download does not use
    # ncbitool (e.g. direct s3 link), indicate that there is no versioning.
    input_fasta_name = os.path.basename(ref_file)
    if ref_file.startswith("s3://"):
        input_fasta_local = fetch_from_s3(
            ref_file, dest_dir, auto_unzip, allow_s3mi=True)
        version_number = VERSION_NONE
    elif ref_file.startswith("ftp://"):
        execute_command("cd %s; wget %s" % (dest_dir, ref_file))
        input_fasta_local = os.path.join(dest_dir, input_fasta_name)
        version_number = VERSION_NONE
    else:
        ncbitool_path = install_ncbitool(ncbitool_dest_dir)
        version_number = get_reference_version_number(ncbitool_path, ref_file)
        input_fasta_local = download_ref_local(
            ncbitool_path, ref_file, version_number, dest_dir)
    return input_fasta_local, version_number


def download_ref_remote_with_tool(ncbitool_path, input_fasta_ncbi_path,
                                  version_number, dest_dir, key_path,
                                  remote_username, server_ip):
    cmd = "cd %s; sudo %s file --download --version-num %s %s" % (
        dest_dir, ncbitool_path, version_number, input_fasta_ncbi_path)
    execute_command(remote_command(cmd, key_path, remote_username, server_ip))
    return os.path.join(dest_dir, os.path.basename(input_fasta_ncbi_path))


def download_ref_remote_with_version_any_type(
        ref_file,
        dest_dir,
        local_ncbitool_dest_dir,
        remote_ncbitool_dest_dir,
        key_path,
        remote_username,
        server_ip,
        sudo=False):
    # Get input reference and version number. If download does not use
    # ncbitool (e.g. direct s3 link), indicate that there is no versioning.
    input_fasta_name = os.path.basename(ref_file)
    if ref_file.startswith("s3://"):
        print("Downloading reference from s3...")
        input_fasta_remote = os.path.join(dest_dir, input_fasta_name)
        cmd = "aws s3 cp --quiet %s %s/" % (ref_file, dest_dir)
        execute_command(
            remote_command(cmd, key_path, remote_username, server_ip))
        version_number = VERSION_NONE
    elif ref_file.startswith("ftp://"):
        print("Downloading reference from ftp...")
        input_fasta_remote = os.path.join(dest_dir, input_fasta_name)
        cmd = "cd %s; sudo wget %s" % (dest_dir, ref_file)
        execute_command(
            remote_command(cmd, key_path, remote_username, server_ip))
        version_number = VERSION_NONE
    else:
        local_ncbitool, remote_ncbitool = install_ncbitool(
            local_ncbitool_dest_dir, remote_ncbitool_dest_dir, key_path,
            remote_username, server_ip, sudo)
        version_number = get_reference_version_number(local_ncbitool, ref_file)
        input_fasta_remote = download_ref_remote_with_tool(
            remote_ncbitool, ref_file, version_number, dest_dir, key_path,
            remote_username, server_ip)
    return input_fasta_remote, version_number


def get_host_index_version_file(star_genome):
    genome_dir = os.path.dirname(star_genome)
    print("Getting host index version file...")
    cmd = "aws s3 ls %s/ | grep version.txt" % genome_dir
    version_file = execute_command_with_output(cmd).rstrip().split(" ")[-1]
    return os.path.join(genome_dir, version_file)


def major_version(version):
    m = re.match("(\d+\.\d+).*", version)
    return m.group(1) if m else None


def upload_version_tracker(source_file, output_name, reference_version_number,
                           output_path_s3, indexing_version):
    version_tracker_file = "%s.version.txt" % output_name
    version_json = {
        "name": output_name,
        "source_file": source_file,
        "source_version": reference_version_number,
        "indexing_version": indexing_version,
        "generation_date": datetime.datetime.now().isoformat()
    }
    with open(version_tracker_file, 'wb') as f:
        json.dump(version_json, f)
    cmd = "aws s3 cp --quiet %s %s/" % (version_tracker_file, output_path_s3)
    execute_command(cmd)


def env_set_if_blank(key, value):
    os.environ[key] = os.environ.get(key, value)


# Notes for lineage and taxonomy ID functions:
#
# A hit will normally have (positive) NCBI taxon IDs at all levels of the
# hierarchy, but we use an artificial negative taxon ID if we have determined
#  that the alignment is not specific at the taxonomy level under
# consideration. This happens when a read's multiple reference matches do not
#  agree on taxon ID at the given level.
#
# For example, a read may match 5 references that all belong to different
# species (e.g. Escherichia albertii, Escherichia vulneris, Escherichia coli,
#  ...), but to the same genus (Escherichia). In this case, we use the taxon
# ID for the genus (Escherichia) at the genus-level, but we populate the
# species-level with an artificial negative ID. The artificial ID is defined
# based on a negative base ( INVALID_CALL_BASE_ID), the taxon level (e.g. 2
# for genus), and the valid parent ID (e.g. genus Escherichia's taxon ID).
#
# See also comments under call_hits_m8().


def cleaned_taxid_lineage(taxid_lineage, hit_taxid_str, hit_level_str):
    """Take the taxon lineage and mark meaningless calls with fake taxids."""
    # This assumption is being made in postprocessing
    assert len(taxid_lineage) == 3
    result = [None, None, None]
    hit_tax_level = int(hit_level_str)
    for tax_level, taxid in enumerate(taxid_lineage, 1):
        if tax_level >= hit_tax_level:
            taxid_str = str(taxid)
        else:
            taxid_str = str(
                tax_level * INVALID_CALL_BASE_ID - int(hit_taxid_str))
        result[tax_level - 1] = taxid_str
    return result


def fill_missing_calls(cleaned_lineage):
    """Replace missing calls with virtual taxids as shown in
    fill_missing_calls_tests. Replaces the negative call IDs with a
    calculated negative value that embeds the next higher positive call in it
    to indicate (1) the call is non-specific (negative/artificial) and (2)
    there is a positive specific call above it on the taxonomy tree.

    Ex: with species_id, genus_id, family_id
    (55, -200, 1534) => (55, -200001534, 1534)
    (-100, 5888, -300) => (-100005888, 5888, -300)
    """
    result = list(cleaned_lineage)
    tax_level = len(cleaned_lineage)
    closest_real_hit_just_above_me = -1

    def blank(taxid_int):
        return 0 > taxid_int > INVALID_CALL_BASE_ID

    while tax_level > 0:
        me = int(cleaned_lineage[tax_level - 1])
        if me >= 0:
            closest_real_hit_just_above_me = me
        elif closest_real_hit_just_above_me >= 0 and blank(me):
            result[tax_level - 1] = str(tax_level * INVALID_CALL_BASE_ID -
                                        closest_real_hit_just_above_me)
        tax_level -= 1
    return result


def fill_missing_calls_tests():
    # -200 => -200001534
    assert fill_missing_calls((55, -200, 1534)) == [55, "-200001534", 1534]
    # -100 => -100005888
    assert fill_missing_calls((-100, 5888, -300)) == ["-100005888", 5888, -300]
    # -100 => -100001534, -200 => -200001534
    assert fill_missing_calls((-100, -200,
                               1534)) == ["-100001534", "-200001534", 1534]
    # no change
    assert fill_missing_calls((55, -200, -300)) == [55, -200, -300]


# a bit ad-hoc
fill_missing_calls_tests()


def validate_taxid_lineage(taxid_lineage, hit_taxid_str, hit_level_str):
    cleaned = cleaned_taxid_lineage(taxid_lineage, hit_taxid_str,
                                    hit_level_str)
    return fill_missing_calls(cleaned)


def mark_job_complete(s3_folder, status):
    done_file = "%s.%s" % (AWS_BATCH_JOB_ID, status)
    execute_command("echo '' | aws s3 cp - %s/%s" % (s3_folder, done_file))
