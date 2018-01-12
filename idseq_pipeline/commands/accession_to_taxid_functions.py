from .common import *

# output location
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3').rstrip('/')

# input reference: path relative to ftp://ftp.ncbi.nlm.nih.gov, but we retrieve it using ncbitool
INPUT = "/pub/taxonomy/accession2taxid/...?"
