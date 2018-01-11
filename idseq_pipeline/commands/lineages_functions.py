import os
from .common import *

# output location
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3').rstrip('/')

def make_lineages():
    execute_command("git clone https://github.com/chanzuckerberg/ncbitax2lin.git")
    execute_command("cd ncbitax2lin; make")
    execute_command("aws s3 cp taxid-lineages.db %s/" % OUTPUT_PATH_S3)
    execute_command("aws s3 cp taxid-lineages.csv.gz %s/" % OUTPUT_PATH_S3)
    execute_command("aws s3 cp names.csv.gz %s/" % OUTPUT_PATH_S3)

