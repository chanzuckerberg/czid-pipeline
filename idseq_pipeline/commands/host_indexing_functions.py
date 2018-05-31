import os
import multiprocessing
from .common import *

MAX_STAR_PART_SIZE = 3252010122

# data directories
# from common import ROOT_DIR
DEST_DIR = ROOT_DIR + '/idseq/indexes'  # generated indexes go here

# arguments from environment variables
INPUT_FASTA_S3 = os.environ.get('INPUT_FASTA_S3')
INPUT_GTF_S3 = os.environ.get('INPUT_GTF_S3')
OUTPUT_PATH_S3 = os.environ.get('OUTPUT_PATH_S3').rstrip('/')
HOST_NAME = os.environ.get('HOST_NAME')

# executables
STAR = "STAR"
BOWTIE2_BUILD = "bowtie2-build"

# output names
STAR_INDEX_OUT = 'STAR_genome.tar'
BOWTIE2_INDEX_OUT = 'bowtie2_genome.tar'


### Functions
def split_fasta(fasta_file, max_fasta_part_size):
    fasta_file_list = []
    part_idx = 0
    current_size = 0
    current_output_file_name = "%s.%d" % (fasta_file, part_idx)
    current_output_file = open(current_output_file_name, 'wb')
    fasta_file_list.append(current_output_file_name)
    with open(fasta_file, 'rb') as input_f:
        current_read = input_f.readline()
        for line in input_f:
            # Check if we have to switch different output fasta file
            if current_size > max_fasta_part_size:
                current_output_file.close()
                part_idx += 1
                current_size = 0
                current_output_file_name = "%s.%d" % (fasta_file, part_idx)
                current_output_file = open(current_output_file_name, 'wb')
                fasta_file_list.append(current_output_file_name)

            if line[0] == '>':  # got a new read
                current_output_file.write(current_read)
                current_size += len(current_read)
                current_read = line
            else:
                current_read += line
        current_output_file.write(current_read)
        current_output_file.close()
    return fasta_file_list


def upload_star_index(result_dir, scratch_dir, star_genome_dir_name):
    # archive
    execute_command("tar cvf %s/%s -C %s %s" %
                    (result_dir, STAR_INDEX_OUT, scratch_dir,
                     star_genome_dir_name))
    # copy to S3
    execute_command("aws s3 cp --quiet %s/%s %s/;" %
                    (result_dir, STAR_INDEX_OUT, OUTPUT_PATH_S3))
    # cleanup
    execute_command("cd %s; rm -rf *" % scratch_dir)


def make_star_index(fasta_file, gtf_file, result_dir, scratch_dir, lazy_run):
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
        print("start making STAR index part %d" % i)
        gtf_command_part = ''
        if i == 0 and gtf_file:
            gtf_command_part = '--sjdbGTFfile %s' % gtf_file

        star_genome_part_dir = "%s/part-%d" % (star_genome_dir_name, i)
        star_command_params = [
            'cd', scratch_dir, ';', 'mkdir -p ', star_genome_part_dir, ';',
            STAR, '--runThreadN',
            str(multiprocessing.cpu_count()), '--runMode', 'genomeGenerate',
            gtf_command_part, '--genomeDir', star_genome_part_dir,
            '--genomeFastaFiles', fasta_file_list[i]
        ]
        execute_command(" ".join(star_command_params))
        print("finished making STAR index part %d " % i)
    # record # parts into parts.txt
    execute_command(" echo %d > %s/%s/parts.txt" %
                    (len(fasta_file_list), scratch_dir, star_genome_dir_name))
    t = MyThread(
        target=upload_star_index,
        args=[result_dir, scratch_dir, star_genome_dir_name])
    t.start()
    return t


def make_bowtie2_index(host_name, fasta_file, result_dir, scratch_dir,
                       lazy_run):
    if lazy_run:
        output = os.path.join(result_dir, BOWTIE2_INDEX_OUT)
        if os.path.isfile(output):
            return 1
    bowtie2_genome_dir_name = BOWTIE2_INDEX_OUT.split('.')[0]
    bowtie2_command_params = [
        'cd', scratch_dir, ';'
        'mkdir', bowtie2_genome_dir_name, ';', 'cd', bowtie2_genome_dir_name,
        ';', BOWTIE2_BUILD, fasta_file, host_name
    ]
    execute_command(" ".join(bowtie2_command_params))
    print("finished making bowtie2 index")
    # archive
    execute_command("tar cvf %s/%s -C %s %s" %
                    (result_dir, BOWTIE2_INDEX_OUT, scratch_dir,
                     bowtie2_genome_dir_name))
    # copy to S3
    execute_command("aws s3 cp --quiet %s/%s %s/;" %
                    (result_dir, BOWTIE2_INDEX_OUT, OUTPUT_PATH_S3))
    # cleanup
    execute_command("cd %s; rm -rf *" % scratch_dir)


def make_indexes(version, lazy_run=False):
    # Set up
    input_fasta_name = os.path.basename(INPUT_FASTA_S3)
    host_name = os.path.splitext(input_fasta_name)[0]
    host_dir = os.path.join(DEST_DIR, host_name)
    fasta_dir = os.path.join(host_dir, 'fastas')
    result_dir = os.path.join(host_dir, 'results')
    scratch_dir_star = os.path.join(host_dir, 'scratch_star')
    scratch_dir_bowtie2 = os.path.join(host_dir, 'scratch_bowtie2')
    execute_command("mkdir -p %s %s %s %s" % (host_dir, fasta_dir, result_dir,
                                              scratch_dir_star))
    execute_command("mkdir -p %s %s %s %s" % (host_dir, fasta_dir, result_dir,
                                              scratch_dir_bowtie2))

    input_gtf_local = None
    print(INPUT_GTF_S3)
    if INPUT_GTF_S3:
        input_gtf_local, _version_number = download_ref_local_with_version_any_type(
            INPUT_GTF_S3, fasta_dir, scratch_dir_star, auto_unzip=True)

    input_fasta_local, version_number = download_ref_local_with_version_any_type(
        INPUT_FASTA_S3, fasta_dir, scratch_dir_star, auto_unzip=True)

    # unzip if necessary --- this is only necessary when the data did not come from S3, and should really
    # not happen like this here --- should be a streaming unzip instead, like in the fetch function in common.py
    if os.path.splitext(input_fasta_local)[1] == ".gz":
        execute_command("gunzip -f %s" % input_fasta_local)
        input_fasta_local = os.path.splitext(input_fasta_local)[0]

    # handle lazy_run
    if lazy_run:
        # Download existing files and see what has been done
        command = "aws s3 cp --quiet %s %s --recursive" % (OUTPUT_PATH_S3,
                                                           result_dir)
        execute_command(command)

    # make STAR index
    star_upload_thread = make_star_index(input_fasta_local, input_gtf_local,
                                         result_dir, scratch_dir_star,
                                         lazy_run)

    # make bowtie2 index
    make_bowtie2_index(host_name, input_fasta_local, result_dir,
                       scratch_dir_bowtie2, lazy_run)

    star_upload_thread.join()
    assert not star_upload_thread.exception

    # upload version tracker file
    if not lazy_run:
        upload_version_tracker(INPUT_FASTA_S3, HOST_NAME, version_number,
                               OUTPUT_PATH_S3, version)
