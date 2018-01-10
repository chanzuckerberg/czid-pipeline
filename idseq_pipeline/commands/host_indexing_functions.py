import os
import subprocess
import multiprocessing
from .common import *

MAX_STAR_PART_SIZE = 3252010122

# data directories
ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + '/idseq/indexes' # generated indexes go here

# arguments from environment variables
INPUT_FASTA_S3 = os.environ.get('INPUT_FASTA_S3')
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3').rstrip('/')

# executables
STAR = "STAR"
BOWTIE2_BUILD = "bowtie2-build"

# output names
STAR_INDEX_OUT = 'STAR_genome.tar.gz'
BOWTIE2_INDEX_OUT = 'bowtie2_genome.tar.gz'

### Functions
def split_fasta(fasta_file, max_fasta_part_size):
    fasta_file_list = []
    part_idx = 0; current_size = 0
    current_output_file_name = "%s.%d" % (fasta_file, part_idx)
    current_output_file = open(current_output_file_name, 'wb')
    fasta_file_list.append(current_output_file_name)
    with open(fasta_file, 'rb') as input_f:
        current_read = input_f.readline()
        for line in input_f:
            # Check if we have to switch different output fasta file
            if current_size > max_fasta_part_size:
                current_output_file.close()
                part_idx += 1; current_size = 0
                current_output_file_name = "%s.%d" % (fasta_file, part_idx)
                current_output_file = open(current_output_file_name, 'wb')
                fasta_file_list.append(current_output_file_name)

            if line[0] == '>': # got a new read
                current_output_file.write(current_read)
                current_size += len(current_read)
                current_read = line
            else:
                current_read += line
        current_output_file.write(current_read)
        current_output_file.close()
    return fasta_file_list


def make_star_index(fasta_file, result_dir, scratch_dir, lazy_run):
    if lazy_run:
        output = os.path.join(result_dir, STAR_INDEX_OUT)
        if os.path.isfile(output):
            return 1
    star_genome_dir_name = STAR_INDEX_OUT.split('.')[0]

    # star genome organization
    # STAR_genome/part-${i}, parts.txt
    fasta_file_list = []
    if os.path.getsize(fasta_file) > MAX_STAR_PART_SIZE:
        fasta_file_list = split_fasta(fasta_file, MAX_STAR_PART_SIZE)
    else:
        fasta_file_list.append(fasta_file)

    for i in range(len(fasta_file_list)):
        print "start making STAR index part %d" % i
        star_genome_part_dir = "%s/part-%d" % (star_genome_dir_name, i)
        star_command_params = ['cd', scratch_dir, ';',
                               'mkdir -p ', star_genome_part_dir, ';',
                               STAR, '--runThreadN', str(multiprocessing.cpu_count()),
                               '--runMode', 'genomeGenerate',
                               '--genomeDir', star_genome_part_dir,
                               '--genomeFastaFiles', fasta_file_list[i]]
        execute_command(" ".join(star_command_params))
        print "finished making STAR index part %d " % i
    # record # parts into parts.txt
    execute_command(" echo %d > %s/%s/parts.txt" % (len(fasta_file_list), scratch_dir, star_genome_dir_name))
    # archive and compress
    execute_command("tar czvf %s/%s -C %s %s" % (result_dir, STAR_INDEX_OUT, scratch_dir, star_genome_dir_name))
    # copy to S3
    execute_command("aws s3 cp %s/%s %s/;" % (result_dir, STAR_INDEX_OUT, OUTPUT_PATH_S3))
    # cleanup
    execute_command("cd %s; rm -rf *" % scratch_dir)

def make_bowtie2_index(host_name, fasta_file, result_dir, scratch_dir, lazy_run):
    if lazy_run:
        output = os.path.join(result_dir, BOWTIE2_INDEX_OUT)
        if os.path.isfile(output):
            return 1
    bowtie2_genome_dir_name = BOWTIE2_INDEX_OUT.split('.')[0]
    bowtie2_command_params = ['cd', scratch_dir, ';'
                              'mkdir', bowtie2_genome_dir_name, ';',
                              'cd', bowtie2_genome_dir_name, ';',
                              BOWTIE2_BUILD, fasta_file, host_name]
    execute_command(" ".join(bowtie2_command_params))
    print "finished making bowtie2 index"
    # archive and compress
    execute_command("tar czvf %s/%s -C %s %s" % (result_dir, BOWTIE2_INDEX_OUT, scratch_dir, bowtie2_genome_dir_name))
    # copy to S3
    execute_command("aws s3 cp %s/%s %s/;" % (result_dir, BOWTIE2_INDEX_OUT, OUTPUT_PATH_S3))
    # cleanup
    execute_command("cd %s; rm -rf *" % scratch_dir)

def make_indexes(lazy_run = True):
    input_fasta_name = os.path.basename(INPUT_FASTA_S3)
    host_name = os.path.splitext(input_fasta_name)[0]
    host_dir = os.path.join(DEST_DIR, host_name)
    fasta_dir = os.path.join(host_dir, 'fastas')
    result_dir = os.path.join(host_dir, 'results')
    scratch_dir = os.path.join(host_dir, 'scratch')
    execute_command("mkdir -p %s %s %s %s" % (host_dir, fasta_dir, result_dir, scratch_dir))

    # Download input
    execute_command("aws s3 cp %s %s/" % (INPUT_FASTA_S3, fasta_dir))
    fasta_file = os.path.join(fasta_dir, input_fasta_name)

    if lazy_run:
       # Download existing files and see what has been done
        command = "aws s3 cp %s %s --recursive" % (OUTPUT_PATH_S3, result_dir)
        execute_command(command)

    # make STAR index
    make_star_index(fasta_file, result_dir, scratch_dir, lazy_run)

    # make bowtie2 index
    make_bowtie2_index(host_name, fasta_file, result_dir, scratch_dir, lazy_run)