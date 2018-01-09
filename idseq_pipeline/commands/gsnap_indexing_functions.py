import os
from .common import *

# address and key of GSNAP machine to use
SERVER_IP = os.environ.get('SERVER_IP')
KEY_S3_PATH = os.environ.get('KEY_S3_PATH')

# output location and name
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3')
OUTPUT_NAME = "nt_k16"

# FASTA reference to index
INPUT_FASTA_S3 = "s3://czbiohub-ncbi-store/blast/db/FASTA/nt.gz"

# location of gmapdb, specified during compilation (see https://github.com/juliangehring/GMAP-GSNAP/blob/master/README)
GMAPDB_PATH = "/home/ubuntu/share"

# location of gmap-gsnap executables
GSNAPL_PATH = "/home/ubuntu/bin"

# working directory on machine to use
WORK_DIR = "/home/ubuntu/gmap_build_workdir"

def make_index():
    # get key
    execute_command("aws s3 cp %s %s/" % (KEY_S3_PATH, WORK_DIR))
    key_path = os.path.join(WORK_DIR, os.path.basename(KEY_S3_PATH))
    remote_username = "ubuntu"
    # get input fasta
    execute_command(remote_command("aws s3 cp %s %s/" % (INPUT_FASTA_S3, WORK_DIR)))
    input_fasta = os.path.join(WORK_DIR, os.path.basename(INPUT_FASTA_S3))
    # make index
    indexing_command = "%s/gmap_build -d %s -k 16 %s" % (GSNAPL_PATH, OUTPUT_NAME, input_fasta)
    execute_command(remote_command(indexing_command, key_path, remote_username, SERVER_IP))
    # upload index
    upload_command = "aws s3 cp %s/%s %s/%s/ --recursive" % (GMAPDB_PATH, OUTPUT_NAME, OUTPUT_PATH_S3, OUTPUT_NAME)
    execute_command(remote_command(upload_command, key_path, remote_username, SERVER_IP))
