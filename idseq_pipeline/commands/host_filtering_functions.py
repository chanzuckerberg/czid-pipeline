import os
import multiprocessing
import re
from .common import * #pylint: disable=wildcard-import

# output files
STAR_OUT1 = 'unmapped.star.1.fq'
STAR_OUT2 = 'unmapped.star.2.fq'
STAR_COUNTS_OUT = 'reads_per_gene.star.tab'
PRICESEQFILTER_OUT1 = 'priceseqfilter.unmapped.star.1.fq'
PRICESEQFILTER_OUT2 = 'priceseqfilter.unmapped.star.2.fq'
FQ2FA_OUT1 = 'priceseqfilter.unmapped.star.1.fasta'
FQ2FA_OUT2 = 'priceseqfilter.unmapped.star.2.fasta'
CDHITDUP_OUT1 = 'cdhitdup.priceseqfilter.unmapped.star.1.fasta'
CDHITDUP_OUT2 = 'cdhitdup.priceseqfilter.unmapped.star.2.fasta'
LZW_OUT1 = 'lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta'
LZW_OUT2 = 'lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta'
BOWTIE2_OUT = 'bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.sam'
EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta'
EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta'
EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT3 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.merged.fasta'
GSNAP_FILTER_SAM = 'gsnap_filter.sam'
EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1 = 'unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta'
EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT2 = 'unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta'
EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT3 = 'unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.merged.fasta'
LOGS_OUT_BASENAME = 'log'
STATS_OUT = 'stats.json'
VERSION_OUT = 'versions.json'
PIPELINE_VERSION_OUT = 'pipeline_version.txt'

# arguments from environment variables
INPUT_BUCKET = os.environ.get('INPUT_BUCKET')
FILE_TYPE = os.environ.get('FILE_TYPE')
OUTPUT_BUCKET = os.environ.get('OUTPUT_BUCKET').rstrip('/')
STAR_GENOME = os.environ.get('STAR_GENOME', 's3://czbiohub-infectious-disease/references/human/STAR_genome.tar')
BOWTIE2_GENOME = os.environ.get('BOWTIE2_GENOME', 's3://czbiohub-infectious-disease/references/human/bowtie2_genome.tar')
GSNAP_GENOME = os.environ.get('GSNAP_GENOME', os.path.join(os.path.dirname(STAR_GENOME), 'hg38_pantro5_k16.tar'))

STAR_BOWTIE_VERSION_FILE_S3 = os.environ.get('STAR_BOWTIE_VERSION_FILE_S3', get_host_index_version_file(STAR_GENOME))
DB_SAMPLE_ID = os.environ['DB_SAMPLE_ID']
AWS_BATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID', 'local')
SAMPLE_S3_OUTPUT_POSTFIX = "/%s" % major_version(PIPELINE_VERSION) if PIPELINE_VERSION else ""
SAMPLE_S3_INPUT_PATH = INPUT_BUCKET.rstrip('/')
SAMPLE_S3_OUTPUT_PATH = OUTPUT_BUCKET+ SAMPLE_S3_OUTPUT_POSTFIX
sample_name = SAMPLE_S3_INPUT_PATH[5:].rstrip('/').replace('/', '-')
SAMPLE_DIR = DEST_DIR + '/' + sample_name
FASTQ_DIR = SAMPLE_DIR + '/fastqs'
RESULT_DIR = SAMPLE_DIR + '/results'
SCRATCH_DIR = SAMPLE_DIR + '/scratch'
DEFAULT_LOGPARAMS = {"sample_s3_output_path": SAMPLE_S3_OUTPUT_PATH}

# versioning

# target outputs by task
TARGET_OUTPUTS_SINGLE = {"run_star": [os.path.join(RESULT_DIR, STAR_OUT1)],
                         "run_priceseqfilter": [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1)],
                         "run_fq2fa": [os.path.join(RESULT_DIR, FQ2FA_OUT1)],
                         "run_cdhitdup": [os.path.join(RESULT_DIR, CDHITDUP_OUT1)],
                         "run_lzw": [os.path.join(RESULT_DIR, LZW_OUT1)],
                         "run_bowtie2": [os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1)],
                         "run_gsnap_filter": [os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1)]}
TARGET_OUTPUTS_PAIRED = {"run_star": [os.path.join(RESULT_DIR, STAR_OUT1),
                                      os.path.join(RESULT_DIR, STAR_OUT2)],
                         "run_priceseqfilter": [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1),
                                                os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2)],
                         "run_fq2fa": [os.path.join(RESULT_DIR, FQ2FA_OUT1),
                                       os.path.join(RESULT_DIR, FQ2FA_OUT2)],
                         "run_cdhitdup": [os.path.join(RESULT_DIR, CDHITDUP_OUT1),
                                          os.path.join(RESULT_DIR, CDHITDUP_OUT2)],
                         "run_lzw": [os.path.join(RESULT_DIR, LZW_OUT1),
                                     os.path.join(RESULT_DIR, LZW_OUT2)],
                         "run_bowtie2": [os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1),
                                         os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2),
                                         os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT3)],
                         "run_gsnap_filter": [os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1),
                                              os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT2),
                                              os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT3)]}

# software packages
STAR = "STAR"
PRICESEQ_FILTER = "PriceSeqFilter"
CDHITDUP = "cd-hit-dup"
BOWTIE2 = "bowtie2"
GSNAPL = "gsnapl"

# pipeline configuration
LZW_FRACTION_CUTOFF = 0.45

# convenience functions
def fq2fa(input_fastq, output_fasta):
    execute_command("sed -n '1~4s/^@/>/p;2~4p' <%s >%s" % (input_fastq, output_fasta))

def lzw_fraction(sequence):
    if sequence == "":
        return 0.0
    sequence = sequence.upper()
    dict_size = 0
    dictionary = {}
    # Initialize dictionary with single char
    for c in sequence:
        dict_size += 1
        dictionary[c] = dict_size
    word = ""
    results = []
    for c in sequence:
        wc = word + c
        if dictionary.get(wc):
            word = wc
        else:
            results.append(dictionary[word])
            dict_size += 1
            dictionary[wc] = dict_size
            word = c
    if word != "":
        results.append(dictionary[word])
    return float(len(results))/len(sequence)

def generate_lzw_filtered_single(fasta_file, output_prefix, cutoff_fraction):
    output_read_1 = open(output_prefix + '.1.fasta', 'wb')
    read_1 = open(fasta_file, 'rb')
    count = 0
    filtered = 0
    while True:
        line_r1_header = read_1.readline()
        line_r1_sequence = read_1.readline()
        if line_r1_header and line_r1_sequence:
            fraction_1 = lzw_fraction(line_r1_sequence.rstrip())
            count += 1
            if fraction_1 > cutoff_fraction:
                output_read_1.write(line_r1_header)
                output_read_1.write(line_r1_sequence)
            else:
                filtered += 1
        else:
            break
    print "LZW filter: total reads: %d, filtered reads: %d, kept ratio: %f" % (count, filtered, 1 - float(filtered)/count)
    output_read_1.close()

def generate_lzw_filtered_paired(fasta_file_1, fasta_file_2, output_prefix, cutoff_fraction):
    output_read_1 = open(output_prefix + '.1.fasta', 'wb')
    output_read_2 = open(output_prefix + '.2.fasta', 'wb')
    read_1 = open(fasta_file_1, 'rb')
    read_2 = open(fasta_file_2, 'rb')
    count = 0
    filtered = 0
    while True:
        line_r1_header = read_1.readline()
        line_r1_sequence = read_1.readline()
        line_r2_header = read_2.readline()
        line_r2_sequence = read_2.readline()
        if line_r1_header and line_r1_sequence and line_r2_header and line_r2_sequence:
            fraction_1 = lzw_fraction(line_r1_sequence.rstrip())
            fraction_2 = lzw_fraction(line_r2_sequence.rstrip())
            count += 1
            if fraction_1 > cutoff_fraction and fraction_2 > cutoff_fraction:
                output_read_1.write(line_r1_header)
                output_read_1.write(line_r1_sequence)
                output_read_2.write(line_r2_header)
                output_read_2.write(line_r2_sequence)
            else:
                filtered += 1
        else:
            break
    print "LZW filter: total reads: %d, filtered reads: %d, kept ratio: %f" % (count, filtered, 1 - float(filtered)/count)
    output_read_1.close()
    output_read_2.close()


def generate_unmapped_singles_from_sam(sam_file, output_prefix):
    "Output a single file containing every unmapped read after bowtie2."
    with open(output_prefix + '.1.fasta', 'wb') as output_read_1:
        with open(sam_file, 'rb') as samf:
            # skip headers
            read1 = samf.readline()
            while read1 and read1[0] == '@':
                read1 = samf.readline()
            while read1:
                parts1 = read1.split("\t")
                if parts1[1] == "4": # read unmapped, see https://broadinstitute.github.io/picard/explain-flags.html
                    output_read_1.write(">%s\n%s\n" %(parts1[0], parts1[9]))  # do NOT append /1 to read id
                read1 = samf.readline()


def generate_unmapped_pairs_from_sam_work(output_read_1, output_read_2, output_merged_read, samf):
    # skip headers
    read1 = samf.readline()
    while read1 and read1[0] == '@':
        read1 = samf.readline()
    read2 = samf.readline()
    while read1 and read2:
        parts1 = read1.split("\t")
        parts2 = read2.split("\t")
        if parts1[1] == "77" and parts2[1] == "141": # both parts unmapped, see https://broadinstitute.github.io/picard/explain-flags.html
            output_read_1.write(">%s\n%s\n" %(parts1[0], parts1[9]))
            output_read_2.write(">%s\n%s\n" %(parts2[0], parts2[9]))
            output_merged_read.write(">%s/1\n%s\n" %(parts1[0], parts1[9]))  # append /1 to read id
            output_merged_read.write(">%s/2\n%s\n" %(parts2[0], parts2[9]))  # append /2 to read id
        read1 = samf.readline()
        read2 = samf.readline()


def generate_unmapped_pairs_from_sam(sam_file, output_prefix):
    """Output 1.fasta and 2.fasta containing the unmapped pairs from bowtie2.
    Also output .merged.fasta multiplxing read ids by appending /1 and /2."""
    with open(output_prefix + '.1.fasta', 'wb') as output_read_1:
        with open(output_prefix + '.2.fasta', 'wb') as output_read_2:
            with open(output_prefix + '.merged.fasta', 'wb') as output_merged_read:
                with open(sam_file, 'rb') as samf:
                    generate_unmapped_pairs_from_sam_work(output_read_1, output_read_2, output_merged_read, samf)


# job functions
def run_star_part(output_dir, genome_dir, fastq_files, count_genes=False):
    execute_command("mkdir -p %s" % output_dir)
    star_command_params = ['cd', output_dir, ';', STAR,
                           '--outFilterMultimapNmax', '99999',
                           '--outFilterScoreMinOverLread', '0.5',
                           '--outFilterMatchNminOverLread', '0.5',
                           '--outReadsUnmapped', 'Fastx',
                           '--outFilterMismatchNmax', '999',
                           '--outSAMmode', 'None',
                           '--clip3pNbases', '0',
                           '--runThreadN', str(multiprocessing.cpu_count()),
                           '--genomeDir', genome_dir,
                           '--readFilesIn', " ".join(fastq_files)]
    if fastq_files[0][-3:] == '.gz':
        star_command_params += ['--readFilesCommand', 'zcat']
    if count_genes and os.path.isfile("%s/sjdbList.fromGTF.out.tab" % genome_dir):
        star_command_params += ['--quantMode', 'GeneCounts']
    execute_command(" ".join(star_command_params), os.path.join(output_dir, "Log.progress.out"))


def uncompressed(s3genome):
    if s3genome.endswith(".gz"):
        return s3genome[:-3]
    if s3genome.endswith(".tgz"):
        return s3genome[:-3] + "tar"
    return s3genome


def fetch_genome_work(s3genome):
    genome_name = os.path.basename(s3genome).rstrip(".gz").rstrip(".tar")
    if genome_name not in ("STAR_genome", "bowtie2_genome", "hg38_pantro5_k16"):
        write_to_log("Oh hello interesting new genome {}".format(genome_name))
    genome_dir = os.path.join(REF_DIR, genome_name)
    if not os.path.exists(genome_dir):
        try:
            #TODO Clean up and fold into fetch_reference
            # Hmm. I understand the desire to create a single function for fetching any reference we might need,
            # but the functionality of these two functions is actually quite different.  One of them handles
            # resiliently situations like the destination existing/not existing, being a folder/not a folder, etc.
            # The other expands an archive into an existing folder, and needs to thoroughly cleanup in case
            # something does not go as planned, because of the prefetching pattern of use.  Also, this function
            # will fetch whichever of the "tar" or "tar.gz" versions of the archive actually exist in s3,
            # even if the other version has been specified.  So the two functions are resilient in different ways.
            install_s3mi()
            tarfile = uncompressed(s3genome)
            try:
                execute_command("s3mi cat {tarfile} | tar xvf - -C {refdir}".format(tarfile=tarfile, refdir=REF_DIR))
                assert os.path.isdir(genome_dir)
            except:
                if tarfile != s3genome:
                    # The uncompressed version doesn't exist.   This is much slower, but no choice.
                    execute_command("rm -rf {}".format(genome_dir))
                    execute_command("s3mi cat {s3genome} | tar xvfz - -C {refdir}".format(s3genome=s3genome, refdir=REF_DIR))
                    assert os.path.isdir(genome_dir)
                else:
                    # Okay, may be s3mi is broken.  We'll try aws cp next.
                    raise
        except:
            try:
                execute_command("rm -rf {}".format(genome_dir))
                execute_command("aws s3 cp --quiet {s3genome} - | tar xvf - -C {refdir}".format(s3genome=s3genome, refdir=REF_DIR))
                assert os.path.isdir(genome_dir)
            except:
                write_to_log("Failed to download index {}, it might not exist.".format(s3genome))
                genome_dir = None
                # Note we do not reraise the exception here, just print it.
                traceback.print_exc()
        if genome_dir:
            write_to_log("downloaded index {}".format(s3genome))
    return genome_dir


def fetch_genome(s3genome, mutex=threading.RLock(), mutexes={}): #pylint: disable=dangerous-default-value
    """Fetch and expand genome archive from s3 into local dir.  Return that local dir.
    If all retries fail, return None.  If already downloaded, return right away.
    If a fetch of the same genome is already in progress on another thread, wait for it.

    Typical use:

            fruitfly_dir = fetch_genome("s3://fruitfly_genome.tar")
            assert fruitfly_dir != None
            # Prefetching optimization:  While doing compute intensive work on fruit flies,
            # start fetching the butterfly genome.
            threading.Thread(target=fetch_genome, args=["s3://butterfly_genome.tar"]).start()
            ... do some compute intensive work on fruit flies ...
            butterfly_dir = fetch_genome("s3://butterfly_genome.tar")
            assert butterfly_dir != None
            threading.Thread(target=fetch_genome, args=["s3://firefly_genome.tar"]).start()
            ... do some compute intensive work on butterflies ...
            firefly_dir = fetch_genome("s3://firefly_genome.tar")
            ...

    Without the prefetching thread, the compute intensive work on butterflies would have
    to wait for the entire butterfly genome to be downloaded.   With prefetching like this,
    the download of the butterfly genome proceeds in parallel with the fruit fly computation,
    and by the time the butterfly genome is needed, it may already have been fully downloaded,
    so the butterfly computation won't have to wait for it.  Similarly, the download of the
    firefly genome proceeds in parallel with the butterfly computation, and the firefly genome
    will be ready by the time it's needed.  The program would still work corerctly if we
    comment out all the prefetching threads, but would take much longer to execute.

    If the different stages of this computation support result caching, so that typically
    the first N stages would be cached from a previous run, and the computation would
    have to resume from stage N+1, this pattern works beautifully by avoiding any unnecessary
    fetching of data that won't be needed for the cached stages, while still fetching the
    data needed for stage N+1.
    """
    with mutex:
        if s3genome not in mutexes:
            mutexes[s3genome] = threading.RLock()
        mx = mutexes[s3genome]
    with mx:
        return fetch_genome_work(s3genome)


def get_read(f):
    # The FASTQ format specifies that each read consists of 4 lines,
    # the first of which begins with @ followed by read ID.
    r, rid = [], None
    line = f.readline()
    if line:
        assert line[0] == "@"
        rid = line.split("\t", 1)[0].strip()
        r.append(line)
        r.append(f.readline())
        r.append(f.readline())
        r.append(f.readline())
    return r, rid


def write_lines(of, lines):
    for l in lines:
        of.write(l)


def handle_outstanding_read(r0, r0id, outstanding_r0, outstanding_r1, of0, of1, mem, max_mem):
    # If read r0 completes an outstanding r1, output the pair (r0, r1).
    # Else r0 becomes outstanding, so in future some r1 may complete it.
    if r0id:
        if r0id in outstanding_r1:
            write_lines(of0, r0)
            write_lines(of1, outstanding_r1.pop(r0id))
            mem -= 1
        else:
            outstanding_r0[r0id] = r0
            mem += 1
            if mem > max_mem:
                max_mem = mem
    return mem, max_mem


def sync_pairs_work(of0, of1, if0, if1):
    # TODO:  Use this as a template for merging fasta?
    outstanding_r0 = {}
    outstanding_r1 = {}
    mem = 0
    max_mem = 0
    while True:
        r0, r0id = get_read(if0)
        r1, r1id = get_read(if1)
        if not r0 and not r1:
            break
        if r0id == r1id:
            # If the input pairs are already synchronized, we take this branch on every iteration.
            write_lines(of0, r0)
            write_lines(of1, r1)
        else:
            mem, max_mem = handle_outstanding_read(r0, r0id, outstanding_r0, outstanding_r1, of0, of1, mem, max_mem)
            mem, max_mem = handle_outstanding_read(r1, r1id, outstanding_r1, outstanding_r0, of1, of0, mem, max_mem)
    return outstanding_r0, outstanding_r1, max_mem


def sync_pairs(fastq_files, max_discrepancies=0):
    """
    The given fastq_files contain the same read IDs but in different order.
    Output the same data in sycnhronized order.  Omit up to max_discrepancies
    if necessary.  If more must be suppressed, raise assertion.
    """
    if len(fastq_files) != 2:
        return fastq_files
    output_filenames = [ifn + ".synchornized_pairs.fq" for ifn in fastq_files]
    with open(fastq_files[0], "rb") as if_0:
        with open(fastq_files[1], "rb") as if_1:
            with open(output_filenames[0], "wb") as of_0:
                with open(output_filenames[1], "wb") as of_1:
                    outstanding_r0, outstanding_r1, max_mem = sync_pairs_work(of_0, of_1, if_0, if_1)
    if max_mem:
        # This will be printed if some pairs were out of order.
        warning_message = "WARNING:  Pair order out of sync in {fqf}.  Synchronized using RAM for {max_mem} pairs.".format(fqf=fastq_files, max_mem=max_mem)
        write_to_log(warning_message)
    discrepancies_count = len(outstanding_r0) + len(outstanding_r1)
    if discrepancies_count:
        warning_message = "WARNING:  Found {dc} broken pairs in {fqf}, e.g., {example}.".format(dc=discrepancies_count, fqf=fastq_files, example=(outstanding_r0 or outstanding_r1).popitem()[0])
        write_to_log(warning_message)
        assert discrepancies_count <= max_discrepancies, warning_message
    return output_filenames


def run_star(fastq_files):
    star_outputs = [STAR_COUNTS_OUT, STAR_OUT1, STAR_OUT2]
    num_fastqs = len(fastq_files)
    gene_count_output = None
    def unmapped_files_in(some_dir):
        return ["%s/Unmapped.out.mate%d" % (some_dir, i+1) for i in range(num_fastqs)]
    genome_dir = fetch_genome(STAR_GENOME)
    assert genome_dir != None
    # If we are here, we are also going to need a bowtie genome later;  start fetching it now
    # This is the absolute PERFECT PLACE for this fetch.  If we are computing from scratch,
    # the download has plenty of time to complete before bowtie needs it.  If we are doing
    # a lazy rerun, this function gets skipped, and we avoid a huge unnecessary download.
    threading.Thread(target=fetch_genome, args=[BOWTIE2_GENOME]).start()
    # Check if parts.txt file exists, if so use the new version of (partitioned indices). Otherwise, stay put
    parts_file = os.path.join(genome_dir, "parts.txt")

    assert os.path.isfile(parts_file)
    with open(parts_file, 'rb') as parts_f:
        num_parts = int(parts_f.read())
    unmapped = fastq_files
    for part_idx in range(num_parts):
        tmp_result_dir = "%s/star-part-%d" % (SCRATCH_DIR, part_idx)
        genome_part = "%s/part-%d" % (genome_dir, part_idx)
        run_star_part(tmp_result_dir, genome_part, unmapped, part_idx == 0)
        unmapped = sync_pairs(unmapped_files_in(tmp_result_dir))
        # run part 0 in gene-counting mode:
        # (a) ERCCs are doped into part 0 and we want their counts
        # (b) if there is only 1 part (e.g. human), the host gene counts also make sense
        if part_idx == 0:
            gene_count_file = os.path.join(tmp_result_dir, "ReadsPerGene.out.tab")
            if os.path.isfile(gene_count_file):
                gene_count_output = gene_count_file

    result_files = [gene_count_output] + unmapped
    for i, f in enumerate(result_files):
        if f != None:
            output_i = os.path.join(RESULT_DIR, star_outputs[i])
            execute_command("mv %s %s;" % (f, output_i))
            execute_command("aws s3 cp --quiet %s %s/;" % (output_i, SAMPLE_S3_OUTPUT_PATH))
    # cleanup
    execute_command("cd %s; rm -rf *" % SCRATCH_DIR)
    write_to_log("finished job")


def run_priceseqfilter(input_fqs):
    # PriceSeqFilter determines input type based on extension.
    # It will throw an exception if output extension doesn't match input.
    correct_file_extension = os.path.splitext(FILE_TYPE)[0]
    input_files = ["%s.%s" % (fq, correct_file_extension) for fq in input_fqs]
    output_files = ["%s_priceseqfilter_output.%s" % (f, correct_file_extension) for f in input_files]
    for fq, f in zip(input_fqs, input_files):
        execute_command("ln %s %s" % (fq, f))
    priceseq_params = [PRICESEQ_FILTER,
                       '-a', '12',
                       '-rnf', '90',
                       '-log', 'c']
    if len(input_fqs) == 2:
        priceseq_params.extend(['-fp', input_files[0], input_files[1],
                                '-op', output_files[0], output_files[1]])
    else:
        priceseq_params.extend(['-f', input_files[0],
                                '-o', output_files[0]])
    if "fastq" in FILE_TYPE:
        priceseq_params.extend(['-rqf', '85', '0.98'])
    execute_command(" ".join(priceseq_params))
    write_to_log("finished job")
    execute_command("mv %s %s" % (output_files[0], os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1)))
    execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, PRICESEQFILTER_OUT1, SAMPLE_S3_OUTPUT_PATH))
    if len(input_fqs) == 2:
        execute_command("mv %s %s" % (output_files[1], os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2)))
        execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, PRICESEQFILTER_OUT2, SAMPLE_S3_OUTPUT_PATH))

def run_fq2fa(input_fqs):
    fq2fa(input_fqs[0], os.path.join(RESULT_DIR, FQ2FA_OUT1))
    if len(input_fqs) == 2:
        fq2fa(input_fqs[1], os.path.join(RESULT_DIR, FQ2FA_OUT2))
    write_to_log("finished job")
    execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, FQ2FA_OUT1, SAMPLE_S3_OUTPUT_PATH))
    if len(input_fqs) == 2:
        execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, FQ2FA_OUT2, SAMPLE_S3_OUTPUT_PATH))

def run_cdhitdup(input_fas):
    cdhitdup_params = [CDHITDUP,
                       '-i', input_fas[0],
                       '-o', RESULT_DIR + '/' + CDHITDUP_OUT1,
                       '-e', '0.05', '-u', '70']
    if len(input_fas) == 2:
        cdhitdup_params.extend(['-i2', input_fas[1],
                                '-o2', RESULT_DIR + '/' + CDHITDUP_OUT2])
    execute_command(" ".join(cdhitdup_params))
    write_to_log("finished job")
    execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, CDHITDUP_OUT1, SAMPLE_S3_OUTPUT_PATH))
    if len(input_fas) == 2:
        execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, CDHITDUP_OUT2, SAMPLE_S3_OUTPUT_PATH))

def run_lzw(input_fas):
    output_prefix = RESULT_DIR + '/' + LZW_OUT1[:-8]
    if len(input_fas) == 2:
        generate_lzw_filtered_paired(input_fas[0], input_fas[1], output_prefix, LZW_FRACTION_CUTOFF)
    else:
        generate_lzw_filtered_single(input_fas[0], output_prefix, LZW_FRACTION_CUTOFF)
    write_to_log("finished job")
    # copy back to aws
    execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, LZW_OUT1, SAMPLE_S3_OUTPUT_PATH))
    if len(input_fas) == 2:
        execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, LZW_OUT2, SAMPLE_S3_OUTPUT_PATH))

def run_bowtie2(input_fas):
    # check if genome downloaded already
    genome_dir = fetch_genome(BOWTIE2_GENOME)
    assert genome_dir != None

    # If we are here, we are also going to need a gsnap genome later;  start fetching it now.
    # This is actually THE PERFECT PLACE to initiate this fetch.  When we are running from
    # scratch, there is plenty of time to download the gsnap genome while bowtie is running.
    # When we are doing a lazy rerun, this function gets skipped, and the fetching of gsnap
    # genome is not initiated.  That's brilliant -- we don't fetch the gsnap genome if
    # we won't be needing it, and lazy reruns are very quick.
    threading.Thread(target=fetch_genome, args=[GSNAP_GENOME]).start()
    # the file structure looks like "bowtie2_genome/GRCh38.primary_assembly.genome.3.bt2"
    # the code below will handle up to "bowtie2_genome/GRCh38.primary_assembly.genome.99.bt2" but not 100
    local_genome_dir_ls = execute_command_with_output("ls {genome_dir}/*.bt2*".format(genome_dir=genome_dir))
    genome_basename = local_genome_dir_ls.split("\n")[0][:-6]
    if genome_basename[-1] == '.':
        genome_basename = genome_basename[:-1]
    bowtie2_params = [BOWTIE2,
                      '-q',
                      '-p', str(multiprocessing.cpu_count()),
                      '-x', genome_basename,
                      '-f',
                      '--very-sensitive-local',
                      '-S', RESULT_DIR + '/' + BOWTIE2_OUT]
    if len(input_fas) == 2:
        bowtie2_params.extend(['-1', input_fas[0], '-2', input_fas[1]])
    else:
        bowtie2_params.extend(['-U', input_fas[0]])
    execute_command(" ".join(bowtie2_params))
    write_to_log("finished alignment")
    # extract out unmapped files from sam
    output_prefix = RESULT_DIR + '/' + EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1[:-8]
    if len(input_fas) == 2:
        generate_unmapped_pairs_from_sam(RESULT_DIR + '/' + BOWTIE2_OUT, output_prefix)
    else:
        generate_unmapped_singles_from_sam(RESULT_DIR + '/' + BOWTIE2_OUT, output_prefix)
    write_to_log("extracted unmapped fragments from SAM file")
    execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, BOWTIE2_OUT, SAMPLE_S3_OUTPUT_PATH))
    execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1, SAMPLE_S3_OUTPUT_PATH))
    if len(input_fas) == 2:
        execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2, SAMPLE_S3_OUTPUT_PATH))
        execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT3, SAMPLE_S3_OUTPUT_PATH))


# Can remove this once the todo below, issue #173, is addressed.
class SkipGsnap(Exception):
    pass


def run_gsnap_filter(input_fas):
    # Unpack the gsnap genome
    genome_dir = fetch_genome(GSNAP_GENOME)
    if genome_dir == None:
        # Apparently if the GSNAP_GENOME file doesn't exist, we are supposed to skip this step.
        # TODO (yunfang):  An independent way to specify whether this step should be executed,
        # so that operational errors don't just silently cause the step to be skipped. See #173.
        raise SkipGsnap()
    gsnap_base_dir = os.path.dirname(genome_dir)
    gsnap_index_name = os.path.basename(genome_dir)

    # Run Gsnap
    gsnap_params = [GSNAPL,
                    '-A sam',
                    '--batch=0',
                    '--use-shared-memory=0',
                    '--gmap-mode=all',
                    '--npaths=1',
                    '--ordered',
                    '-t 32',
                    '--max-mismatches=40',
                    '-D', gsnap_base_dir,
                    '-d', gsnap_index_name,
                    '-o', RESULT_DIR + '/' + GSNAP_FILTER_SAM]
    gsnap_params += input_fas
    execute_command(" ".join(gsnap_params))
    write_to_log("finished alignment")
    # extract out unmapped files from sam
    output_prefix = RESULT_DIR + '/' + EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1[:-8]
    if len(input_fas) == 2:
        generate_unmapped_pairs_from_sam(RESULT_DIR + '/' + GSNAP_FILTER_SAM, output_prefix)
    else:
        generate_unmapped_singles_from_sam(RESULT_DIR + '/' + GSNAP_FILTER_SAM, output_prefix)
    write_to_log("extracted unmapped fragments from SAM file for gsnap output")
    execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, GSNAP_FILTER_SAM, SAMPLE_S3_OUTPUT_PATH))
    execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1, SAMPLE_S3_OUTPUT_PATH))
    if len(input_fas) == 2:
        execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT2, SAMPLE_S3_OUTPUT_PATH))
        execute_command("aws s3 cp --quiet %s/%s %s/;" % (RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT3, SAMPLE_S3_OUTPUT_PATH))


def run_host_filtering(fastq_files, initial_file_type_for_log, lazy_run, stats):
    number_of_input_files = len(fastq_files)
    target_outputs = TARGET_OUTPUTS_PAIRED if number_of_input_files == 2 else TARGET_OUTPUTS_SINGLE

    # run STAR
    logparams = return_merged_dict(
        DEFAULT_LOGPARAMS,
        {"title": "STAR",
         "version_file_s3": STAR_BOWTIE_VERSION_FILE_S3,
         "output_version_file": os.path.join(RESULT_DIR, VERSION_OUT)})
    run_and_log_s3(logparams, target_outputs["run_star"], lazy_run, SAMPLE_S3_OUTPUT_PATH, run_star, fastq_files)
    stats.count_reads("run_star",
                      before_filename=fastq_files[0],
                      before_filetype=initial_file_type_for_log,
                      after_filename=os.path.join(RESULT_DIR, STAR_OUT1),
                      after_filetype=initial_file_type_for_log)

    # run priceseqfilter
    logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "PriceSeqFilter"})
    if number_of_input_files == 2:
        input_files = [os.path.join(RESULT_DIR, STAR_OUT1), os.path.join(RESULT_DIR, STAR_OUT2)]
    else:
        input_files = [os.path.join(RESULT_DIR, STAR_OUT1)]
    run_and_log_s3(logparams, target_outputs["run_priceseqfilter"], lazy_run, SAMPLE_S3_OUTPUT_PATH, run_priceseqfilter, input_files)
    stats.count_reads("run_priceseqfilter",
                      before_filename=os.path.join(RESULT_DIR, STAR_OUT1),
                      before_filetype=initial_file_type_for_log,
                      after_filename=os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1),
                      after_filetype=initial_file_type_for_log)

    # run fastq to fasta
    if "fastq" in FILE_TYPE:
        logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "FASTQ to FASTA"})
        if number_of_input_files == 2:
            input_files = [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1), os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2)]
            next_inputs = [os.path.join(RESULT_DIR, FQ2FA_OUT1), os.path.join(RESULT_DIR, FQ2FA_OUT2)]
        else:
            input_files = [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1)]
            next_inputs = [os.path.join(RESULT_DIR, FQ2FA_OUT1)]
        run_and_log_s3(logparams, target_outputs["run_fq2fa"], lazy_run, SAMPLE_S3_OUTPUT_PATH, run_fq2fa, input_files)
    else:
        if number_of_input_files == 2:
            next_inputs = [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1), os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2)]
        else:
            next_inputs = [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1)]

    # run cdhitdup
    logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "CD-HIT-DUP"})
    run_and_log_s3(logparams, target_outputs["run_cdhitdup"], lazy_run, SAMPLE_S3_OUTPUT_PATH, run_cdhitdup, next_inputs)
    stats.count_reads("run_cdhitdup",
                      before_filename=next_inputs[0],
                      before_filetype="fasta_paired",
                      after_filename=os.path.join(RESULT_DIR, CDHITDUP_OUT1),
                      after_filetype="fasta_paired")

    # run lzw filter
    logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "LZW filter"})
    if number_of_input_files == 2:
        input_files = [os.path.join(RESULT_DIR, CDHITDUP_OUT1), os.path.join(RESULT_DIR, CDHITDUP_OUT2)]
    else:
        input_files = [os.path.join(RESULT_DIR, CDHITDUP_OUT1)]
    run_and_log_s3(logparams, target_outputs["run_lzw"], lazy_run, SAMPLE_S3_OUTPUT_PATH, run_lzw, input_files)
    stats.count_reads("run_lzw",
                      before_filename=os.path.join(RESULT_DIR, CDHITDUP_OUT1),
                      before_filetype="fasta_paired",
                      after_filename=os.path.join(RESULT_DIR, LZW_OUT1),
                      after_filetype="fasta_paired")

    # run bowtie
    logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "bowtie2"})
    if number_of_input_files == 2:
        input_files = [os.path.join(RESULT_DIR, LZW_OUT1), os.path.join(RESULT_DIR, LZW_OUT2)]
    else:
        input_files = [os.path.join(RESULT_DIR, LZW_OUT1)]
    run_and_log_s3(logparams, target_outputs["run_bowtie2"], lazy_run, SAMPLE_S3_OUTPUT_PATH, run_bowtie2, input_files)
    stats.count_reads("run_bowtie2",
                      before_filename=os.path.join(RESULT_DIR, LZW_OUT1),
                      before_filetype="fasta_paired",
                      after_filename=os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1),
                      after_filetype="fasta_paired")

    # run gsnap against host genomes (only available for human as of 5/1/2018)
    try:
        if number_of_input_files == 2:
            input_files = [os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1), os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2)]
        else:
            input_files = [os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1)]
        logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "run_gsnap_filter"})
        run_and_log_s3(logparams, target_outputs["run_gsnap_filter"], lazy_run, SAMPLE_S3_OUTPUT_PATH, run_gsnap_filter, input_files)
        stats.count_reads("run_gsnap_filter",
                          before_filename=os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1),
                          before_filetype="fasta_paired",
                          after_filename=os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1),
                          after_filetype="fasta_paired")
    except SkipGsnap:
        pass

    # finalize the remaing reads
    stats.set_remaining_reads()


def upload_pipeline_version_file():
    execute_command("echo %s > %s" % (PIPELINE_VERSION, PIPELINE_VERSION_OUT))
    execute_command("aws s3 cp %s %s/" % (PIPELINE_VERSION_OUT, OUTPUT_BUCKET))


def run_stage1(lazy_run=True):
    execute_command("mkdir -p %s %s %s %s" % (SAMPLE_DIR, FASTQ_DIR, RESULT_DIR, SCRATCH_DIR))
    execute_command("mkdir -p %s " % REF_DIR)

    # configure logger
    log_file = "%s/%s.%s.txt" % (RESULT_DIR, LOGS_OUT_BASENAME, AWS_BATCH_JOB_ID)
    configure_logger(log_file)

    # Download fastqs
    command = "aws s3 ls %s/ | grep '\\.%s$'" % (SAMPLE_S3_INPUT_PATH, FILE_TYPE)
    output = execute_command_with_output(command).rstrip().split("\n")
    for line in output:
        m = re.match(".*?([^ ]*." + re.escape(FILE_TYPE) + ")", line)
        if m:
            execute_command("aws s3 cp --quiet %s/%s %s/" % (SAMPLE_S3_INPUT_PATH, m.group(1), FASTQ_DIR))
        else:
            print "%s doesn't match %s" % (line, FILE_TYPE)
    fastq_files = execute_command_with_output("ls %s/*.%s" % (FASTQ_DIR, FILE_TYPE)).rstrip().split("\n")

    # Identify input files and characteristics
    if len(fastq_files) not in [1, 2]:
        write_to_log("Number of input files was neither 1 nor 2. Aborting computation.")
        return # only support either 1 file or 2 (paired) files

    # Record total number of input reads
    initial_file_type_for_log = "fastq" if "fastq" in FILE_TYPE else "fasta"
    if len(fastq_files) == 2:
        initial_file_type_for_log += "_paired"
    stats = StatsFile(STATS_OUT, RESULT_DIR, None, SAMPLE_S3_OUTPUT_PATH)
    # TODO: When we start flowing prefiltered jobs, that is, jobs which specify total reads in the input,
    # do not overwrite it by doing this.
    stats.data.append({'total_reads': count_reads(fastq_files[0], initial_file_type_for_log)})

    # run host filtering
    run_host_filtering(fastq_files, initial_file_type_for_log, lazy_run, stats)

    # This lets the webapp know the stage has completed.
    stats.save_to_s3()

    write_to_log("Host filtering complete")
    upload_log_file(SAMPLE_S3_OUTPUT_PATH)
