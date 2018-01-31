import shelve
import subprocess
import os
import gzip
import sys
from .common import *

# arguments from environment variables
INPUT_FASTA_S3 = os.environ.get('INPUT_FASTA_S3')
ACCESSION2TAXID_DB_S3_PATH =  os.environ.get('ACCESSION2TAXID_DB_S3_PATH')

# data directories
ROOT_DIR = '/mnt/idseq'
REF_DIR = ROOT_DIR + '/ref'

# Functions
def make_blacklist():
    execute_command("mkdir -p %s " % REF_DIR)

    # Download accession-to-taxid shelf:
    accession2taxid_gz = os.path.basename(ACCESSION2TAXID_DB_S3_PATH)
    accession2taxid_path = os.path.join(REF_DIR, accession2taxid_gz[:-3])
    if not os.path.isfile(accession2taxid_path):
        execute_command("aws s3 cp %s %s/" % (ACCESSION2TAXID_DB_S3_PATH, REF_DIR))
        execute_command("cd %s; gunzip -f %s" % (REF_DIR, accession2taxid_gz))

    # Generate blacklist:
    accession2taxid_dict = shelve.open(accession2taxid_path)
    execute_command("mkdir -p %s/data; aws s3 cp %s %s/data/" % (ROOT_DIR, INPUT_FASTA_S3, ROOT_DIR))
    input_gz_file = "%s/data/%s" % (ROOT_DIR, os.path.basename(INPUT_FASTA_S3))
    with gzip.open(input_gz_file, 'rb') as f:
        for line in f:
            if line[0] == '>':
                accession_id = line.split('|')[3]
                accession_main = accession_id.split('.')[0]
                taxon_id = accession2taxid_dict.get(accession_main, '-1')
                print ",".join([accession_id, taxon_id])
