import os
from .common import *

# address and key of RAPSearch machine to use
SERVER_IP = os.environ.get('SERVER_IP')
KEY_S3_PATH = os.environ.get('KEY_S3_PATH')

# output location and name
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3').rstrip('/')
OUTPUT_NAME = os.environ.get('OUTPUT_NAME')

# path to gzipped FASTA reference to index, relative to ftp://ftp.ncbi.nlm.nih.gov
INPUT_FASTA_S3 = os.environ.get('INPUT_FASTA_S3')

# location of index builder executable
PRE_RAPSEARCH_PATH = "/usr/local/bin/prerapsearch"

# S3 location of ncbitool executable
NCBITOOL_S3_PATH = "s3://idseq-database/ncbitool"

# working directories
WORK_DIR = "/home/ec2-user/data/pre_rapsearch_workdir"  # on RAPSearch machine
REMOTE_USERNAME = "ec2-user"
KEY_PATH = None
LOCAL_WORK_DIR = "idseq_pipeline_temp"  # locally


def get_key():
    global KEY_PATH
    execute_command("aws s3 cp --quiet %s %s/" % (KEY_S3_PATH, LOCAL_WORK_DIR))
    KEY_PATH = os.path.join(LOCAL_WORK_DIR, os.path.basename(KEY_S3_PATH))
    execute_command("chmod 400 %s" % KEY_PATH)


def make_index(version):
    # Set up
    execute_command("mkdir -p %s" % LOCAL_WORK_DIR)
    get_key()
    execute_command(
        remote_command("sudo mkdir -p %s" % WORK_DIR, KEY_PATH,
                       REMOTE_USERNAME, SERVER_IP))
    local_ncbitool, remote_ncbitool = install_ncbitool(
        LOCAL_WORK_DIR, WORK_DIR, KEY_PATH, REMOTE_USERNAME, SERVER_IP)

    # download reference and unzip
    input_fasta_zipped, version_number = download_reference_on_remote_with_version_any_source_type(
        INPUT_FASTA_S3, WORK_DIR, LOCAL_WORK_DIR, WORK_DIR, KEY_PATH,
        REMOTE_USERNAME, SERVER_IP)
    input_fasta_unzipped = input_fasta_zipped[:-3]  # gunzip removes .gz
    command = "sudo gunzip -f %s" % input_fasta_zipped
    execute_command(
        remote_command(command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

    # make index
    remote_output_dir = os.path.join(WORK_DIR, OUTPUT_NAME)
    execute_command(
        remote_command("sudo mkdir -p %s" % remote_output_dir, KEY_PATH,
                       REMOTE_USERNAME, SERVER_IP))
    indexing_command = "cd %s; sudo %s -d %s -n %s" % (remote_output_dir,
                                                       PRE_RAPSEARCH_PATH,
                                                       input_fasta_unzipped,
                                                       OUTPUT_NAME)
    execute_command(
        remote_command(indexing_command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

    # upload index
    upload_command = "aws s3 cp --quiet %s %s/%s/ --recursive" % (
        remote_output_dir, OUTPUT_PATH_S3, OUTPUT_NAME)
    execute_command(
        remote_command(upload_command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

    # upload input fasta unzipped
    upload_command = "aws s3 cp --quiet %s %s/" % (input_fasta_unzipped,
                                                   OUTPUT_PATH_S3)
    execute_command(
        remote_command(upload_command, KEY_PATH, REMOTE_USERNAME, SERVER_IP))

    # upload version tracker file
    upload_version_tracker(INPUT_FASTA_S3, OUTPUT_NAME, version_number,
                           OUTPUT_PATH_S3, version)
