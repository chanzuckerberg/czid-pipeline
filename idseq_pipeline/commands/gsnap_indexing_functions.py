import os
from .common import *

# address and key of GSNAP machine to use
SERVER_IP = os.environ.get('SERVER_IP')
KEY_S3_PATH = os.environ.get('KEY_S3_PATH')

# output location and name
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3')
OUTPUT_NAME = os.environ.get('OUTPUT_NAME')

# path to gzipped FASTA reference to index, relative to ftp://ftp.ncbi.nlm.nih.gov
INPUT_FASTA_S3 = os.environ.get('INPUT_FASTA_S3')

# location of gmapdb, specified during compilation (see https://github.com/juliangehring/GMAP-GSNAP/blob/master/README)
GMAPDB_PATH = "/home/ubuntu/share"

# location of gmap-gsnap executables
GSNAPL_PATH = "/home/ubuntu/bin"

# S3 location of ncbitool executable
NCBITOOL_S3_PATH = "s3://czbiohub-infectious-disease/ncbitool"

# working directories
WORK_DIR = "/home/ubuntu/gmap_build_workdir" # on GSNAP machine
REMOTE_USERNAME = "ubuntu"
KEY_PATH = None
LOCAL_WORK_DIR = "idseq_pipeline_temp" # locally

def get_key():
    global KEY_PATH
    execute_command("aws s3 cp %s %s/" % (KEY_S3_PATH, LOCAL_WORK_DIR))
    KEY_PATH = os.path.join(LOCAL_WORK_DIR, os.path.basename(KEY_S3_PATH))
    execute_command("chmod 400 %s" % KEY_PATH)

def install_ncbitool():
    # install locally
    execute_command("aws s3 cp %s %s/" % (NCBITOOL_S3_PATH, LOCAL_WORK_DIR))
    execute_command("chmod u+x %s/ncbitool" % LOCAL_WORK_DIR)
    # install on remote
    command = "aws s3 cp %s %s/; " % (NCBITOOL_S3_PATH, WORK_DIR)
    command += "chmod u+x %s/ncbitool" % WORK_DIR
    execute_command(remote_command(command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

def get_reference_version_number():
    command = "%s/ncbitool file history %s" % (LOCAL_WORK_DIR, INPUT_FASTA_S3)
    output = execute_command_with_output(command).split("File History: ")[1]
    version_history = json.loads(output)
    version_numbers = [entry["Version"] for entry in version_history]
    return max(version_numbers) # use latest available version

def download_reference_on_remote(version_number):
    command = "%s/ncbitool file --download --version-num %s %s %s/" % (WORK_DIR, version_number, INPUT_FASTA_S3, WORK_DIR)
    execute_command(remote_command(command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))
    return os.path.join(WORK_DIR, os.path.basename(INPUT_FASTA_S3))

def make_index():
    # Set up
    execute_command("mkdir -p %s" % LOCAL_WORK_DIR)
    get_key()
    execute_command(remote_command("mkdir -p %s" % WORK_DIR, KEY_PATH, REMOTE_USERNAME, SERVER_IP))
    install_ncbitool()

    # Get latest version of reference and unzip
    version_number = get_reference_version_number()
    input_fasta_zipped = download_reference_on_remote(version_number)
    input_fasta_unzipped = input_fasta_zipped[:-3]
    command = "rm -f %s; gunzip %s" % (input_fasta_unzipped, input_fasta_zipped)
    execute_command(remote_command(command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

    # make index
    indexing_command = "sudo %s/gmap_build -d %s -k 16 %s" % (GSNAPL_PATH, OUTPUT_NAME, input_fasta_unzipped)
    execute_command(remote_command(indexing_command, key_path, remote_username, SERVER_IP))

    # upload index
    upload_command = "aws s3 cp %s/%s %s/%s/ --recursive" % (GMAPDB_PATH, OUTPUT_NAME, OUTPUT_PATH_S3, OUTPUT_NAME)
    execute_command(remote_command(upload_command, key_path, remote_username, SERVER_IP))
