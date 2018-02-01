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
WORK_DIR = "/home/ubuntu/share/gmap_build_workdir" # on GSNAP machine
REMOTE_USERNAME = "ubuntu"
KEY_PATH = None
LOCAL_WORK_DIR = "idseq_pipeline_temp" # locally

def get_key():
    global KEY_PATH
    execute_command("aws s3 cp --quiet %s %s/" % (KEY_S3_PATH, LOCAL_WORK_DIR))
    KEY_PATH = os.path.join(LOCAL_WORK_DIR, os.path.basename(KEY_S3_PATH))
    execute_command("chmod 400 %s" % KEY_PATH)

def make_index(version):
    # Set up
    execute_command("mkdir -p %s" % LOCAL_WORK_DIR)
    get_key()
    execute_command(remote_command("sudo mkdir -p %s" % WORK_DIR, KEY_PATH, REMOTE_USERNAME, SERVER_IP))
    local_ncbitool, remote_ncbitool = install_ncbitool(LOCAL_WORK_DIR, WORK_DIR, KEY_PATH, REMOTE_USERNAME, SERVER_IP, True)

    # get latest version number of desired reference file
    version_number = get_reference_version_number(local_ncbitool, INPUT_FASTA_S3)

    # download reference and unzip
    input_fasta_zipped = download_reference_on_remote(remote_ncbitool, INPUT_FASTA_S3, version_number, WORK_DIR,
                                                      KEY_PATH, REMOTE_USERNAME, SERVER_IP)
    input_fasta_unzipped = input_fasta_zipped[:-3]
    command = "sudo gunzip -f %s" % input_fasta_zipped
    execute_command(remote_command(command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

    # make index
    indexing_command = "sudo %s/gmap_build -d %s -k 16 %s" % (GSNAPL_PATH, OUTPUT_NAME, input_fasta_unzipped)
    execute_command(remote_command(indexing_command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

    # upload index
    upload_command = "aws s3 cp --quiet %s/%s %s/%s/ --recursive" % (GMAPDB_PATH, OUTPUT_NAME, OUTPUT_PATH_S3, OUTPUT_NAME)
    execute_command(remote_command(upload_command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

    # upload version tracker file
    upload_version_tracker(INPUT_FASTA_S3, OUTPUT_NAME, version_number, OUTPUT_PATH_S3, version)
