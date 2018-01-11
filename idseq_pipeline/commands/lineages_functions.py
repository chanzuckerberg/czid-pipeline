import os
from .common import *

# output location
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3').rstrip('/')

# input reference: path relative to ftp://ftp.ncbi.nlm.nih.gov, but we retrieve it using ncbitool
INPUT = "/pub/taxonomy/taxdump.tar.gz"

def make_lineages():
    # Install ncbitax2lin
    execute_command("git clone https://github.com/chanzuckerberg/ncbitax2lin.git")
    execute_command("cd ncbitax2lin")
    work_dir = os.getcwd()

    # Get input reference and version number.
    ncbitool_path = install_ncbitool(work_dir)
    version_number = get_reference_version_number(ncbitool_path, INPUT)
    input_fasta_local = download_reference_locally(ncbitool_path, INPUT, version_number, work_dir)
    execute_command("cd %s; rm -rfv taxdump; mkdir -p taxdump/taxdump; tar zxf %s -C ./taxdump/taxdump" % (work_dir, os.path.basename(INPUT)))

    # Make lineage files
    execute_command("cd %s; make" % work_dir)
    execute_command("aws s3 cp taxid-lineages.db %s/" % OUTPUT_PATH_S3)
    execute_command("aws s3 cp taxid-lineages.csv.gz %s/" % OUTPUT_PATH_S3)
    execute_command("aws s3 cp names.csv.gz %s/" % OUTPUT_PATH_S3)

    # Make deuterostome list
