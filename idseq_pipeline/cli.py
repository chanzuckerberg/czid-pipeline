"""
idseq_pipeline

Usage:
  idseq_pipeline host_filtering
  idseq_pipeline non_host_alignment
  idseq_pipeline postprocess
  idseq_pipeline host_indexing
  idseq_pipeline gsnap_indexing
  idseq_pipeline blacklist
  idseq_pipeline -h | --help
  idseq_pipeline --version

Options:
  -h --help                         Show this screen.
  --version                         Show version.

Examples:

  # "Arguments" are set via environment variables:
  FILE_TYPE="fastq.gz" DB_SAMPLE_ID="77" INPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/fastqs" OUTPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/results" FASTQ_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/fastqs" idseq_pipeline host_filtering

  FILE_TYPE="fastq.gz" DB_SAMPLE_ID="77" INPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/results" OUTPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/results" FASTQ_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/fastqs" idseq_pipeline non_host_alignment

  INPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/results" OUTPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/postprocess" idseq_pipeline postprocess

  INPUT_FASTA_S3=s3://czbiohub-infectious-disease/references/mosquitos/mosquito_genomes2.fa OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references/mosquitos idseq_pipeline host_indexing

  INPUT_FASTA_S3=/blast/db/FASTA/nt.gz SERVER_IP=34.211.67.166 KEY_S3_PATH=s3://czbiohub-infectious-disease/idseq-production.pem OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references OUTPUT_NAME=nt_k16 idseq_pipeline gsnap_indexing

  INPUT_FASTA_S3='s3://czbiohub-ncbi-store/blast/db/FASTA/vector.gz' ACCESSION2TAXID_DB_S3_PATH='s3://czbiohub-infectious-disease/references/accession2taxid.db.gz' idseq_pipeline blacklist

Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/chanzuckerberg/idseq-web
"""


from inspect import getmembers, isclass

from docopt import docopt

from . import __version__ as VERSION


def main():
    """Main CLI entrypoint."""
    import idseq_pipeline.commands
    options = docopt(__doc__, version=VERSION)

    # Match user-entered command
    for (k, v) in options.items(): 
        if hasattr(idseq_pipeline.commands, k) and v:
            module = getattr(idseq_pipeline.commands, k)
            idseq_pipeline.commands = getmembers(module, isclass)
            command = [command[1] for command in idseq_pipeline.commands if command[0] != 'Base'][0]
            command = command(options)
            command.run()
