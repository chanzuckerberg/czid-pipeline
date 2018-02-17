"""
idseq_pipeline

Usage:
  idseq_pipeline host_filtering
  idseq_pipeline non_host_alignment
  idseq_pipeline postprocess
  idseq_pipeline host_indexing
  idseq_pipeline gsnap_indexing
  idseq_pipeline rapsearch_indexing
  idseq_pipeline blacklist
  idseq_pipeline lineages
  idseq_pipeline curate_accession2taxid --mapping_files <mapping_file1,mapping_file2,etc> --nr_file <nr_file> --nt_file <nt_file> --output_s3_folder <output_s3_folder> [--previous_mapping <previous_mapping>]
  idseq_pipeline curate_accessionid2seq --input_s3_path <nt_or_nr_file> --output_s3_path <output_db>
  idseq_pipeline accessionid2seq --s3_db_path <ccession2seq_db> --input_annotated_fasta_path <fasta> --ouput_path <output_json_path>
  idseq_pipeline -h | --help
  idseq_pipeline --version

Options:
  -h --help                                             Show this screen.
  --version                                             Show version.
  --mapping_files <mapping_file1,mapping_file2,etc>     Accession2taxid files from NCBI, comma-separated [default: /pub/taxonomy/accession2taxid/nucl_est.accession2taxid.gz,/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz,/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz,/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz,/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz,/pub/taxonomy/accession2taxid/prot.accession2taxid.gz]
  --nr_file <nr_file>                                   NR fasta file [default: /blast/db/FASTA/nr.gz]
  --nt_file <nt_file>                                   NT fasta file [default: /blast/db/FASTA/nt.gz]
  --output_s3_folder <output_s3_folder>                 S3 destination for curated accesstion2taxid output file [default: s3://czbiohub-infectious-disease/references]
  --previous_mapping <previous_mapping>                 S3 path to the previous version of accesstion2taxid.db.gz, to be compared to the new version

Examples:

  # "Arguments" are set via environment variables:
  FILE_TYPE="fastq.gz" DB_SAMPLE_ID="77" INPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/fastqs" OUTPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/results" FASTQ_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/fastqs" idseq_pipeline host_filtering

  FILE_TYPE="fastq.gz" DB_SAMPLE_ID="77" INPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/results" OUTPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/results" FASTQ_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/fastqs" idseq_pipeline non_host_alignment

  INPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/results" OUTPUT_BUCKET="s3://czbiohub-idseq-samples-development/samples/3/77/postprocess" idseq_pipeline postprocess

  INPUT_FASTA_S3=s3://czbiohub-infectious-disease/references/mosquitos/mosquito_genomes2.fa OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references/mosquitos idseq_pipeline host_indexing

  INPUT_FASTA_S3=ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references/human idseq_pipeline host_indexing

  INPUT_FASTA_S3=/blast/db/FASTA/nt.gz SERVER_IP=34.211.67.166 KEY_S3_PATH=s3://czbiohub-infectious-disease/idseq-production.pem OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references OUTPUT_NAME=nt_k16 idseq_pipeline gsnap_indexing

  INPUT_FASTA_S3=/blast/db/FASTA/nr.gz SERVER_IP=54.191.193.210 KEY_S3_PATH=s3://czbiohub-infectious-disease/idseq-alpha.pem OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references OUTPUT_NAME=nr_rapsearch idseq_pipeline rapsearch_indexing

  OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references idseq_pipeline lineages

  INPUT_FASTA_S3='s3://czbiohub-ncbi-store/blast/db/FASTA/vector.gz' ACCESSION2TAXID_DB_S3_PATH='s3://czbiohub-infectious-disease/references/accession2taxid.db.gz' idseq_pipeline blacklist

Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/chanzuckerberg/idseq-web
"""


from inspect import getmembers, isclass

from docopt import docopt

from version import __version__

def main():
    """Main CLI entrypoint."""
    import commands as idseq_pipeline_commands
    options = docopt(__doc__, version=__version__)

    # Match user-entered command
    for (k, v) in options.items():
        if hasattr(idseq_pipeline_commands, k) and v:
            module = getattr(idseq_pipeline_commands, k)
            pipeline_commands = getmembers(module, isclass)
            command = [command[1] for command in pipeline_commands if command[0] != 'Base'][0]
            command = command(options)
            command.version = __version__
            command.run()
