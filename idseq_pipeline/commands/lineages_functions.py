import os
from .common import *

# output location
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3').rstrip('/')

# input reference: path relative to ftp://ftp.ncbi.nlm.nih.gov, but we retrieve it using ncbitool
INPUT = "/pub/taxonomy/taxdump.tar.gz"

def make_lineages():
    # Install ncbitax2lin
    execute_command("rm -rf ncbitax2lin; git clone https://github.com/chanzuckerberg/ncbitax2lin.git")
    work_dir = os.path.join(os.getcwd(), "ncbitax2lin")

    # Get input reference and version number
    ncbitool_path = install_ncbitool(work_dir)
    version_number = get_reference_version_number(ncbitool_path, INPUT)
    input_fasta_local = download_reference_locally(ncbitool_path, INPUT, version_number, work_dir)
    command = "cd %s; rm -rfv taxdump; mkdir -p taxdump/taxdump; " % work_dir
    command += "tar zxf %s -C ./taxdump/taxdump" % os.path.basename(INPUT) # this structure is needed for ncbitax2lin 'make'
    execute_command(command)

    # Make and upload lineage files
    execute_command("cd %s; make" % work_dir)
    execute_command("aws s3 cp %s/taxid-lineages.db %s/" % (work_dir, OUTPUT_PATH_S3))
    execute_command("aws s3 cp %s/taxid-lineages.csv.gz %s/" % (work_dir, OUTPUT_PATH_S3))
    execute_command("aws s3 cp %s/names.csv.gz %s/" % (work_dir, OUTPUT_PATH_S3))

    # Make and upload deuterostome list
    input_filename = "lineages.csv"
    output_filename = "deuterostome_taxids.txt"
    command = "cd %s; gunzip %s.gz; grep 'Deuterostomia' %s | cut -f1 -d',' > %s; " % (work_dir, input_filename, input_filename, output_filename)
    command += "aws s3 cp %s %s/" % (output_filename, OUTPUT_PATH_S3)
    execute_command(command)

    # upload version:
    upload_version_tracker('lineage_and_deuterostome', version_number, OUTPUT_PATH_S3)
