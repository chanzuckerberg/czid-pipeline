#!/bin/bash

## Install idseq-pipeline
pip install git+https://github.com/chanzuckerberg/idseq-pipeline.git

## Archive path
DATE=`date '+%Y-%m-%d_%H-%M-%S'`
ARCHIVE_FOLDER=s3://czbiohub-infectious-disease/references/archive/$DATE

## Make GSNAP index
# (a) archive old index
aws s3 cp s3://czbiohub-infectious-disease/references/nt_k16.tar ${ARCHIVE_FOLDER}/
aws s3 cp s3://czbiohub-infectious-disease/references/nt_k16.version.txt ${ARCHIVE_FOLDER}/
# (b) make and upload new index
INPUT_FASTA_S3=/blast/db/FASTA/nt.gz SERVER_IP=34.211.67.166 KEY_S3_PATH=s3://czbiohub-infectious-disease/idseq-production.pem OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references OUTPUT_NAME=nt_k16 idseq_pipeline gsnap_indexing

## Make RAPSearch2 index
# (a) archive old index
aws s3 cp s3://czbiohub-infectious-disease/references/nr_rapsearch/nr_rapsearch ${ARCHIVE_FOLDER}/
aws s3 cp s3://czbiohub-infectious-disease/references/nr_rapsearch/nr_rapsearch.info ${ARCHIVE_FOLDER}/
# (b) make and upload new index
INPUT_FASTA_S3=/blast/db/FASTA/nr.gz SERVER_IP=54.191.193.210 KEY_S3_PATH=s3://czbiohub-infectious-disease/idseq-alpha.pem OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references OUTPUT_NAME=nr_rapsearch idseq_pipeline rapsearch_indexing

  OUTPUT_PATH_S3=s3://czbiohub-infectious-disease/references idseq_pipeline lineages

  idseq_pipeline curate_accession2taxid --mapping_files <mapping_file1,mapping_file2,etc> --nr_file <nr_file> --nt_file <nt_file> --output_s3_folder <output_s3_folder> [--previous_mapping <previous_mapping>]

