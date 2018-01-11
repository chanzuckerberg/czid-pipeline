import time
import datetime
import threading
import sys
import subprocess
import logging
import json
import gzip
import os
from .. import __version__

NCBITOOL_S3_PATH = "s3://czbiohub-infectious-disease/ncbitool" # S3 location of ncbitool executable

STATS = []
LOGGER = None

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
        super(CommandTracker, self).__init__(update_period, self.print_update)
        with CommandTracker.lock:
            self.id = CommandTracker.count
            CommandTracker.count += 1

    def print_update(self, t_elapsed):
        print "Command %d still running after %3.1f seconds." % (self.id, t_elapsed)
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


def execute_command_with_output(command, progress_file=None):
    with CommandTracker() as ct:
        print "Command {}: {}".format(ct.id, command)
        with ProgressFile(progress_file):
            return subprocess.check_output(command, shell=True)


def execute_command_realtime_stdout(command, progress_file=None):
    with CommandTracker() as ct:
        print "Command {}: {}".format(ct.id, command)
        with ProgressFile(progress_file):
            subprocess.check_call(command, shell=True)


def execute_command(command, progress_file=None):
    execute_command_realtime_stdout(command, progress_file)

def remote_command(base_command, key_path, remote_username, instance_ip):
    return 'ssh -o "StrictHostKeyChecking no" -i %s %s@%s "%s"' % (key_path, remote_username, instance_ip, base_command)

class TimeFilter(logging.Filter):
    def filter(self, record):
        try:
            last = self.last
        except AttributeError:
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

def count_reads(file_name, file_type):
    count = 0
    if file_name[-3:] == '.gz':
        f = gzip.open(file_name)
    else:
        f = open(file_name)
    for line in f:
        if file_type == "fastq_paired":
            count += 2./4
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
    global LOGGER
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

def load_existing_stats(stats_file):
    global STATS
    if os.path.isfile(stats_file):
        with open(stats_file) as f:
            STATS = json.load(f)

def get_total_reads_from_stats():
    for item in STATS:
        if "total_reads" in item:
            return item["total_reads"]
    # check run star if total reads not available
    for item in STATS:
        if item.get("task") == 'run_star':
            return item["reads_before"]
    # no entry
    return 0.1

def get_remaining_reads_from_stats():
    return (item for item in STATS if item.get("task") == "run_gsnapl_remotely").next().get("reads_before")

def run_and_log(logparams, target_outputs, lazy_run, func_name, *args):
    global LOGGER
    LOGGER = logging.getLogger()
    LOGGER.info("========== %s ==========" % logparams.get("title"))
    # copy log file -- start
    LOGGER.handlers[0].flush()
    execute_command("aws s3 cp %s %s/;" % (LOGGER.handlers[0].baseFilename, logparams["sample_s3_output_path"]))
    # produce the output
    if lazy_run and all(os.path.isfile(output) for output in target_outputs):
        LOGGER.info("output exists, lazy run")
    else:
        func_name(*args)
        LOGGER.info("uploaded output")
    # copy log file -- after work is done
    execute_command("aws s3 cp %s %s/;" % (LOGGER.handlers[0].baseFilename, logparams["sample_s3_output_path"]))
    # count records
    required_params = ["before_file_name", "before_file_type", "after_file_name", "after_file_type"]
    if logparams.get("count_reads") and all(param in logparams for param in required_params):
        records_before = count_reads(logparams["before_file_name"], logparams["before_file_type"])
        records_after = count_reads(logparams["after_file_name"], logparams["after_file_type"])
        STATS.append({'task': func_name.__name__, 'reads_before': records_before, 'reads_after': records_after})
    # copy log file -- end
    LOGGER.handlers[0].flush()
    execute_command("aws s3 cp %s %s/;" % (LOGGER.handlers[0].baseFilename, logparams["sample_s3_output_path"]))
    # write stats
    stats_path = logparams.get("stats_file")
    if stats_path and os.path.isfile(stats_path):
        with open(stats_path, 'wb') as f:
            json.dump(STATS, f)
        execute_command("aws s3 cp %s %s/;" % (stats_path, logparams["sample_s3_output_path"]))

def write_to_log(message):
    LOGGER.info(message)

def unbuffer_stdout():
    # Unbuffer stdout and redirect stderr into stdout. This helps observe logged events in realtime.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    os.dup2(sys.stdout.fileno(), sys.stderr.fileno())

def upload_commit_sha():
    sha_file = os.environ.get('COMMIT_SHA_FILE')
    s3_destination = os.environ.get('OUTPUT_BUCKET')
    if sha_file is None or s3_destination is None:
        return
    sha_file_parts = os.path.splitext(os.path.basename(sha_file))
    aws_batch_job_id = os.environ.get('AWS_BATCH_JOB_ID', 'local')
    sha_file_new_name = "%s_job-%s%s" % (sha_file_parts[0], aws_batch_job_id, sha_file_parts[1])
    execute_command("aws s3 cp %s %s/%s;" % (sha_file, s3_destination.rstrip('/'), sha_file_new_name))

def install_ncbitool_locally(local_work_dir):
    execute_command("aws s3 cp %s %s/" % (NCBITOOL_S3_PATH, local_work_dir))
    execute_command("chmod u+x %s/ncbitool" % local_work_dir)
    return "%s/ncbitool" % local_work_dir

def install_ncbitool(local_work_dir, remote_work_dir=None, key_path=None, remote_username=None, server_ip=None):
    local_result = install_ncbitool_locally(local_work_dir)
    if remote_work_dir is None:
        return local_result
    command = "sudo aws s3 cp %s %s/; " % (NCBITOOL_S3_PATH, remote_work_dir)
    command += "sudo chmod u+x %s/ncbitool" % remote_work_dir
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

def download_reference_on_remote(ncbitool_path, input_fasta_ncbi_path, version_number, destination_dir,
                                 key_path, remote_username, server_ip):
    command = "cd %s; sudo %s file --download --version-num %s %s" % (destination_dir, ncbitool_path, version_number, input_fasta_ncbi_path)
    execute_command(remote_command(command, key_path, remote_username, server_ip))
    return os.path.join(destination_dir, os.path.basename(input_fasta_ncbi_path))

def upload_version_tracker(output_name, reference_version_number, output_path_s3):
    version_tracker_file = "%s.version.txt" % output_name
    execute_command("echo 'reference: %d' > %s" % (reference_version_number, version_tracker_file))
    execute_command("echo 'indexing: %s' > %s" % (__version__, version_tracker_file))
    execute_command("aws s3 cp %s %s/" % (version_tracker_file, output_path_s3))
