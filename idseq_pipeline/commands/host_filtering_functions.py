import os
import multiprocessing
import re
from .common import *  #pylint: disable=wildcard-import

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
STATS_IN = 'total_reads.json'
STATS_OUT = 'stats.json'
VERSION_OUT = 'versions.json'
PIPELINE_VERSION_OUT = 'pipeline_version.txt'
MAX_INPUT_READS = 75 * 1000 * 1000
MAX_GNSAP_FILTER_READS = 10 * 1000 * 1000
INPUT_TRUNCATED_FILE = 'input_truncated.txt'

# arguments from environment variables
INPUT_BUCKET = os.environ.get('INPUT_BUCKET')
FILE_TYPE = os.environ.get('FILE_TYPE')
OUTPUT_BUCKET = os.environ.get('OUTPUT_BUCKET').rstrip('/')
STAR_GENOME = os.environ.get(
    'STAR_GENOME',
    's3://czbiohub-infectious-disease/references/human/STAR_genome.tar')
BOWTIE2_GENOME = os.environ.get(
    'BOWTIE2_GENOME',
    's3://czbiohub-infectious-disease/references/human/bowtie2_genome.tar')
GSNAP_GENOME = os.environ.get('GSNAP_GENOME',
                              os.path.join(
                                  os.path.dirname(STAR_GENOME),
                                  'hg38_pantro5_k16.tar'))

STAR_BOWTIE_VERSION_FILE_S3 = os.environ.get(
    'STAR_BOWTIE_VERSION_FILE_S3', get_host_index_version_file(STAR_GENOME))
DB_SAMPLE_ID = os.environ['DB_SAMPLE_ID']
AWS_BATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID', 'local')
SAMPLE_S3_OUTPUT_POSTFIX = "/%s" % major_version(
    PIPELINE_VERSION) if PIPELINE_VERSION else ""
SAMPLE_S3_INPUT_PATH = INPUT_BUCKET.rstrip('/')
SAMPLE_S3_OUTPUT_PATH = OUTPUT_BUCKET + SAMPLE_S3_OUTPUT_POSTFIX
sample_name = SAMPLE_S3_INPUT_PATH[5:].rstrip('/').replace('/', '-')
SAMPLE_DIR = DEST_DIR + '/' + sample_name
FASTQ_DIR = SAMPLE_DIR + '/fastqs'
RESULT_DIR = SAMPLE_DIR + '/results'
SCRATCH_DIR = SAMPLE_DIR + '/scratch'
DEFAULT_LOG_PARAMS = {"sample_s3_output_path": SAMPLE_S3_OUTPUT_PATH}

# versioning

# target outputs by task
TARGET_OUTPUTS_SINGLE = {
    "run_star": [os.path.join(RESULT_DIR, STAR_OUT1)],
    "run_priceseqfilter": [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1)],
    "run_fq2fa": [os.path.join(RESULT_DIR, FQ2FA_OUT1)],
    "run_cdhitdup": [os.path.join(RESULT_DIR, CDHITDUP_OUT1)],
    "run_lzw": [os.path.join(RESULT_DIR, LZW_OUT1)],
    "run_bowtie2":
    [os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1)],
    "run_gsnap_filter":
    [os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1)]
}
TARGET_OUTPUTS_PAIRED = {
    "run_star":
    [os.path.join(RESULT_DIR, STAR_OUT1),
     os.path.join(RESULT_DIR, STAR_OUT2)],
    "run_priceseqfilter": [
        os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1),
        os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2)
    ],
    "run_fq2fa": [
        os.path.join(RESULT_DIR, FQ2FA_OUT1),
        os.path.join(RESULT_DIR, FQ2FA_OUT2)
    ],
    "run_cdhitdup": [
        os.path.join(RESULT_DIR, CDHITDUP_OUT1),
        os.path.join(RESULT_DIR, CDHITDUP_OUT2)
    ],
    "run_lzw":
    [os.path.join(RESULT_DIR, LZW_OUT1),
     os.path.join(RESULT_DIR, LZW_OUT2)],
    "run_bowtie2": [
        os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1),
        os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2),
        os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT3)
    ],
    "run_gsnap_filter": [
        os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1),
        os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT2),
        os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT3)
    ]
}

# software packages
STAR = "STAR"
PRICESEQ_FILTER = "PriceSeqFilter"
CDHITDUP = "cd-hit-dup"
BOWTIE2 = "bowtie2"
GSNAPL = "gsnapl"

# pipeline configuration
LZW_FRACTION_CUTOFFS = [0.45, 0.42]


# Convenience functions
def fq2fa(input_fastq, output_fasta):
    cmd = "sed -n '1~4s/^@/>/p;2~4p' <%s >%s" % (input_fastq, output_fasta)
    execute_command(cmd)


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
    return float(len(results)) / len(sequence)


def generate_lzw_filtered_single(fasta_file, out_prefix, cutoff_fractions):
    out_read_1 = open(out_prefix + '.1.fasta', 'wb')
    for cutoff_frac in cutoff_fractions:
        read_1 = open(fasta_file, 'rb')
        count = 0
        filtered = 0
        while True:
            line_r1_header = read_1.readline()
            line_r1_sequence = read_1.readline()
            if line_r1_header and line_r1_sequence:
                fraction_1 = lzw_fraction(line_r1_sequence.rstrip())
                count += 1
                if fraction_1 > cutoff_frac:
                    out_read_1.write(line_r1_header)
                    out_read_1.write(line_r1_sequence)
                else:
                    filtered += 1
            else:
                break
        msg = "LZW filter: cutoff_frac: %f, total reads: %d, filtered reads: %d, " \
              "kept ratio: %f" % (cutoff_frac, count, filtered, 1 - float(filtered) / count)
        print(msg)
        if count != filtered:
            break
    out_read_1.close()


def generate_lzw_filtered_paired(fasta_file_1, fasta_file_2, out_prefix,
                                 cutoff_fractions):
    out_read_1 = open(out_prefix + '.1.fasta', 'wb')
    out_read_2 = open(out_prefix + '.2.fasta', 'wb')
    for cutoff_frac in cutoff_fractions:
        read_1 = open(fasta_file_1, 'rb')
        read_2 = open(fasta_file_2, 'rb')
        count = 0
        filtered = 0
        while True:
            line_r1_header = read_1.readline()
            line_r1_seq = read_1.readline()
            line_r2_header = read_2.readline()
            line_r2_seq = read_2.readline()
            if line_r1_header and line_r1_seq and line_r2_header and line_r2_seq:
                fraction_1 = lzw_fraction(line_r1_seq.rstrip())
                fraction_2 = lzw_fraction(line_r2_seq.rstrip())
                count += 1
                if fraction_1 > cutoff_frac and fraction_2 > cutoff_frac:
                    out_read_1.write(line_r1_header)
                    out_read_1.write(line_r1_seq)
                    out_read_2.write(line_r2_header)
                    out_read_2.write(line_r2_seq)
                else:
                    filtered += 1
            else:
                break
        msg = "LZW filter: cutoff_frac: %f, total reads: %d, filtered reads: %d, " \
              "kept ratio: %f" % (cutoff_frac, count, filtered, 1 - float(filtered) / count)
        print(msg)
        if count != filtered:
            break
    out_read_1.close()
    out_read_2.close()


def generate_unmapped_singles_from_sam(sam_file, output_prefix):
    """Output a single file containing every unmapped read after bowtie2.

    SAM file alignments:
    - See: https://en.wikipedia.org/wiki/SAM_(file_format)
         - https://broadinstitute.github.io/picard/explain-flags.html
    - part[0] = query template name
    - part[1] = bitwise flag
    - part[9] = segment sequence
    """
    with open(output_prefix + '.1.fasta', 'wb') as output_read:
        with open(sam_file, 'rb') as sam_f:
            # Skip headers
            read = sam_f.readline()
            while read and read[0] == '@':
                read = sam_f.readline()
            while read:
                part = read.split("\t")
                # Read unmapped
                if part[1] == "4":
                    # Do NOT append /1 to read id
                    output_read.write(">%s\n%s\n" % (part[0], part[9]))
                read = sam_f.readline()


def generate_unmapped_pairs_from_sam_work(out_read_1, out_read_2,
                                          out_merged_read, sam_f):
    """SAM file alignments:
    - See: https://en.wikipedia.org/wiki/SAM_(file_format)
         - https://broadinstitute.github.io/picard/explain-flags.html
    - part[0] = query template name
    - part[1] = bitwise flag
    - part[9] = segment sequence
    """
    # Skip headers
    read1 = sam_f.readline()
    while read1 and read1[0] == '@':
        read1 = sam_f.readline()
    read2 = sam_f.readline()

    while read1 and read2:
        part1 = read1.split("\t")
        part2 = read2.split("\t")
        if part1[1] == "77" and part2[1] == "141":  # Both parts unmapped
            out_read_1.write(">%s\n%s\n" % (part1[0], part1[9]))
            out_read_2.write(">%s\n%s\n" % (part2[0], part2[9]))
            # Append /1 to read id
            out_merged_read.write(">%s/1\n%s\n" % (part1[0], part1[9]))
            # Append /2 to read id
            out_merged_read.write(">%s/2\n%s\n" % (part2[0], part2[9]))
        read1 = sam_f.readline()
        read2 = sam_f.readline()


def generate_unmapped_pairs_from_sam(sam_file, out_prefix):
    """Output 1.fasta and 2.fasta containing the unmapped pairs from bowtie2.
    Also output .merged.fasta and multiplex read ids by appending /1 and /2.
    """
    with open(out_prefix + '.1.fasta', 'wb') as out_read_1:
        with open(out_prefix + '.2.fasta', 'wb') as out_read_2:
            with open(out_prefix + '.merged.fasta', 'wb') as out_merged_read:
                with open(sam_file, 'rb') as sam_f:
                    generate_unmapped_pairs_from_sam_work(
                        out_read_1, out_read_2, out_merged_read, sam_f)


def max_input_lines(input_file):
    """Returning number of lines corresponding to MAX_INPUT_READS based on file
    type.
    """
    if "fasta" in input_file:
        return MAX_INPUT_READS * 2
    # Assume it's FASTQ
    return MAX_INPUT_READS * 4


# Job functions


def run_star_part(output_dir, genome_dir, fastq_files, count_genes=False):
    execute_command("mkdir -p %s" % output_dir)
    star_command_params = [
        'cd', output_dir, ';', STAR, '--outFilterMultimapNmax', '99999',
        '--outFilterScoreMinOverLread', '0.5', '--outFilterMatchNminOverLread',
        '0.5', '--outReadsUnmapped', 'Fastx', '--outFilterMismatchNmax', '999',
        '--outSAMmode', 'None', '--clip3pNbases', '0', '--runThreadN',
        str(multiprocessing.cpu_count()), '--genomeDir', genome_dir,
        '--readFilesIn', " ".join(fastq_files)
    ]
    if fastq_files[0][-3:] == '.gz':
        # Create a custom decompressor which does "zcat $input_file | head -
        # ${max_lines}"
        cmd = "echo 'zcat ${2} | head -${1}' > %s/gzhead; " % genome_dir
        execute_command(cmd)
        max_lines = max_input_lines(fastq_files[0])
        star_command_params += [
            '--readFilesCommand',
            '"sh %s/gzhead %d"' % (genome_dir, max_lines)
        ]
    path = "%s/sjdbList.fromGTF.out.tab" % genome_dir
    if count_genes and os.path.isfile(path):
        star_command_params += ['--quantMode', 'GeneCounts']
    cmd = " ".join(star_command_params), os.path.join(output_dir,
                                                      "Log.progress.out")
    execute_command(cmd)


def uncompressed(s3genome):
    if s3genome.endswith(".gz"):
        return s3genome[:-3]
    if s3genome.endswith(".tgz"):
        return s3genome[:-3] + "tar"
    return s3genome


def fetch_genome_work(s3genome, strict):
    genome_name = os.path.basename(s3genome).rstrip(".gz").rstrip(".tar")
    if genome_name not in ("STAR_genome", "bowtie2_genome",
                           "hg38_pantro5_k16"):
        write_to_log("Oh hello interesting new genome {}".format(genome_name))
    genome_dir = os.path.join(REF_DIR, genome_name)

    if not os.path.exists(genome_dir):
        try:
            install_s3mi()
            tarfile = uncompressed(s3genome)
            try:
                print("Trying to download compressed genome...")
                cmd = "s3mi cat {tarfile} | tar xvf - -C {refdir}".format(
                    tarfile=tarfile, refdir=REF_DIR)
                execute_command(cmd)
                assert os.path.isdir(genome_dir)
            except:
                if tarfile != s3genome:
                    print("Uncompressed version doesn't exist. Downloading "
                          "compressed version...")
                    # The uncompressed version doesn't exist. This is much
                    # slower, but no choice.
                    execute_command("rm -rf {}".format(genome_dir))
                    cmd = "s3mi cat {s3genome} | tar xvfz - -C {refdir}".format(
                        s3genome=s3genome, refdir=REF_DIR)
                    execute_command(cmd)
                    assert os.path.isdir(genome_dir)
                else:
                    # Okay, may be s3mi is broken.  We'll try aws cp next.
                    print("Error in downloading with s3mi. Trying aws cp...")
                    raise
        except:
            try:
                execute_command("rm -rf {}".format(genome_dir))
                cmd = "aws s3 cp --quiet {s3genome} - | tar xvf - -C {refdir}".format(
                    s3genome=s3genome, refdir=REF_DIR)
                execute_command(cmd)
                assert os.path.isdir(genome_dir)
            except:
                msg = "Failed to download index {}, it might not exist.".format(
                    s3genome)
                write_to_log(msg)
                if strict:
                    raise
                genome_dir = None
                # Note we do not reraise the exception here, just print it.
                traceback.print_exc()
        if genome_dir:
            write_to_log("successfully downloaded index {}".format(s3genome))
    return genome_dir


def fetch_genome(s3genome, strict=True, mutex=threading.RLock(), mutexes={}):  #pylint: disable=dangerous-default-value
    """Fetch and expand genome archive from s3 into local dir. Return that local
    dir. If already downloaded, return right away. If a fetch of the same genome
    is already in progress on another thread, wait for it to complete or fail;
    and if it failed, try again. If all tries fail, raise an exception (strict)
    or return None (not strict).

    Typical use:

        fruitfly_dir = fetch_genome("s3://fruitfly_genome.tar")
        # Prefetching optimization:  While doing compute intensive work on fruit flies,
        # start fetching the butterfly genome.
        threading.Thread(target=fetch_genome, args=["s3://butterfly_genome.tar"]).start()
        ... do some compute intensive work on fruit flies ...
        butterfly_dir = fetch_genome("s3://butterfly_genome.tar")
        threading.Thread(target=fetch_genome, args=["s3://firefly_genome.tar"]).start()
        ... do some compute intensive work on butterflies ...
        firefly_dir = fetch_genome("s3://firefly_genome.tar")
        ...

    Without the pre-fetching thread, the compute intensive work on butterflies
    would have to wait for the entire butterfly genome to be downloaded. With
    pre-fetching like this, the download of the butterfly genome proceeds in
    parallel with the fruit fly computation, and by the time the butterfly
    genome is needed, it may already have been fully downloaded, so the
    butterfly computation won't have to wait for it. Similarly, the download
    of the firefly genome proceeds in parallel with the butterfly computation,
    and the firefly genome will be ready by the time it's needed. The program
    would still work correctly if we comment out all the pre-fetching threads,
    but would take much longer to execute.

    It may be tempting to initiate all fetching parallel at the start of the
    program, but that's undesirable for two reasons:

      1) Fetching data with s3mi fully utilizes the I/O bandwidth of moderately
    sized instances, so fetching multiple streams in parallel will just slow them
    down and delay the moment we can begin computing.

      2) If the different computation stages support result caching, so that
      typically the first N stages would be cached from a previous run, and the
      computation would resume from stage N+1, this pattern beautifully avoids
      any unnecessary fetching of data that won't be needed for the cached
      stages, while still fetching the data needed for stage N+1. If, instead,
      we were to initiate pre-fetching at the beginning of the program, we would
      have to carefully ensure we only prefetch data that will in fact be
      needed, by replicating some of the caching logic.
    """
    with mutex:
        if s3genome not in mutexes:
            mutexes[s3genome] = threading.RLock()
        mx = mutexes[s3genome]
    with mx:
        return fetch_genome_work(s3genome, strict)


def get_read(f):
    # The FASTQ format specifies that each read consists of 4 lines,
    # the first of which begins with @ followed by read ID.
    read, rid = [], None
    line = f.readline()
    if line:
        assert line[0] == "@"
        rid = line.split("\t", 1)[0].strip()
        read.append(line)
        read.append(f.readline())
        read.append(f.readline())
        read.append(f.readline())
    return read, rid


def write_lines(of, lines):
    for l in lines:
        of.write(l)


def handle_outstanding_read(r0, r0id, outstanding_r0, outstanding_r1, of0, of1,
                            mem, max_mem):
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
            # If the input pairs are already synchronized, we take this branch
            # on every iteration.
            write_lines(of0, r0)
            write_lines(of1, r1)
        else:
            mem, max_mem = handle_outstanding_read(r0, r0id, outstanding_r0,
                                                   outstanding_r1, of0, of1,
                                                   mem, max_mem)
            mem, max_mem = handle_outstanding_read(r1, r1id, outstanding_r1,
                                                   outstanding_r0, of1, of0,
                                                   mem, max_mem)
    return outstanding_r0, outstanding_r1, max_mem


def sync_pairs(fastq_files, max_discrepancies=0):
    """The given fastq_files contain the same read IDs but in different order.
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
                    outstanding_r0, outstanding_r1, max_mem = sync_pairs_work(
                        of_0, of_1, if_0, if_1)
    if max_mem:
        # This will be printed if some pairs were out of order.
        warning_message = "WARNING:  Pair order out of sync in {fqf}.  Synchronized using RAM for {max_mem} pairs.".format(
            fqf=fastq_files, max_mem=max_mem)
        write_to_log(warning_message)
    discrepancies_count = len(outstanding_r0) + len(outstanding_r1)
    if discrepancies_count:
        warning_message = "WARNING:  Found {dc} broken pairs in {fqf}, e.g., {example}.".format(
            dc=discrepancies_count,
            fqf=fastq_files,
            example=(outstanding_r0 or outstanding_r1).popitem()[0])
        write_to_log(warning_message)
        assert discrepancies_count <= max_discrepancies, warning_message
    return output_filenames


def extract_total_counts_from_star_output(result_dir, num_fastqs,
                                          total_counts_from_star):
    ''' Grab the total reads from the Log.final.out file '''
    log_file = os.path.join(result_dir, "Log.final.out")
    total_reads = execute_command_with_output(
        "grep 'Number of input reads' %s" % log_file).split("\t")[1]
    total_reads = int(total_reads)
    if total_reads == MAX_INPUT_READS:  # if it's exactly the same. it must have been truncated
        total_counts_from_star['truncated'] = 1
    total_counts_from_star['total_reads'] = total_reads * num_fastqs


def run_star(fastq_files, uploader_start, total_counts_from_star):
    ''' Run STAR to filter out host '''
    star_outputs = [STAR_COUNTS_OUT, STAR_OUT1, STAR_OUT2]
    num_fastqs = len(fastq_files)
    gene_count_output = None

    def unmapped_files_in(some_dir):
        return [
            "%s/Unmapped.out.mate%d" % (some_dir, i + 1)
            for i in range(num_fastqs)
        ]

    genome_dir = fetch_genome(STAR_GENOME)
    assert genome_dir is not None
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
        # (c) at part 0, we can also extract out total input reads and if the total_counts is exactly the
        #     same as MAX_INPUT_READS then we know the input file is truncated.
        if part_idx == 0:
            gene_count_file = os.path.join(tmp_result_dir,
                                           "ReadsPerGene.out.tab")
            extract_total_counts_from_star_output(tmp_result_dir, num_fastqs,
                                                  total_counts_from_star)
            if os.path.isfile(gene_count_file):
                gene_count_output = gene_count_file

    result_files = [gene_count_output] + unmapped
    for i, f in enumerate(result_files):
        if f is not None:
            output_i = os.path.join(RESULT_DIR, star_outputs[i])
            execute_command("mv %s %s;" % (f, output_i))
            uploader_start(output_i, SAMPLE_S3_OUTPUT_PATH + "/")
    # cleanup
    execute_command("cd %s; rm -rf *" % SCRATCH_DIR)
    write_to_log("finished job")


def run_priceseqfilter(input_fqs, uploader_start):
    # PriceSeqFilter determines input type based on extension.
    # It will throw an exception if output extension doesn't match input.
    correct_file_extension = os.path.splitext(FILE_TYPE)[0]
    input_files = ["%s.%s" % (fq, correct_file_extension) for fq in input_fqs]
    output_files = [
        "%s_priceseqfilter_output.%s" % (f, correct_file_extension)
        for f in input_files
    ]
    for fq, f in zip(input_fqs, input_files):
        execute_command("ln %s %s" % (fq, f))
    priceseq_params = [PRICESEQ_FILTER, '-a', '12', '-rnf', '90', '-log', 'c']
    if len(input_fqs) == 2:
        priceseq_params.extend([
            '-fp', input_files[0], input_files[1], '-op', output_files[0],
            output_files[1]
        ])
    else:
        priceseq_params.extend(['-f', input_files[0], '-o', output_files[0]])
    if "fastq" in FILE_TYPE:
        priceseq_params.extend(['-rqf', '85', '0.98'])
    execute_command(" ".join(priceseq_params))
    write_to_log("finished job")
    execute_command("mv %s %s" %
                    (output_files[0],
                     os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1)))
    uploader_start(
        os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1),
        SAMPLE_S3_OUTPUT_PATH + "/")
    if len(input_fqs) == 2:
        execute_command("mv %s %s" %
                        (output_files[1],
                         os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2)))
        uploader_start(
            os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2),
            SAMPLE_S3_OUTPUT_PATH + "/")


def run_fq2fa(input_fqs, uploader_start):
    fq2fa(input_fqs[0], os.path.join(RESULT_DIR, FQ2FA_OUT1))
    if len(input_fqs) == 2:
        fq2fa(input_fqs[1], os.path.join(RESULT_DIR, FQ2FA_OUT2))
    write_to_log("finished job")
    uploader_start(
        os.path.join(RESULT_DIR, FQ2FA_OUT1), SAMPLE_S3_OUTPUT_PATH + "/")
    if len(input_fqs) == 2:
        uploader_start(
            os.path.join(RESULT_DIR, FQ2FA_OUT2), SAMPLE_S3_OUTPUT_PATH + "/")


def run_cdhitdup(input_fas, uploader_start):
    cdhitdup_params = [
        CDHITDUP, '-i', input_fas[0], '-o', RESULT_DIR + '/' + CDHITDUP_OUT1,
        '-e', '0.05', '-u', '70'
    ]
    if len(input_fas) == 2:
        cdhitdup_params.extend(
            ['-i2', input_fas[1], '-o2', RESULT_DIR + '/' + CDHITDUP_OUT2])
    execute_command(" ".join(cdhitdup_params))
    uploader_start(
        os.path.join(RESULT_DIR, CDHITDUP_OUT1), SAMPLE_S3_OUTPUT_PATH + "/")
    if len(input_fas) == 2:
        uploader_start(
            os.path.join(RESULT_DIR, CDHITDUP_OUT2),
            SAMPLE_S3_OUTPUT_PATH + "/")
    write_to_log("finished job")


def run_lzw(input_fas, uploader_start):
    output_prefix = RESULT_DIR + '/' + LZW_OUT1[:-8]
    if len(input_fas) == 2:
        generate_lzw_filtered_paired(input_fas[0], input_fas[1], output_prefix,
                                     LZW_FRACTION_CUTOFFS)
    else:
        generate_lzw_filtered_single(input_fas[0], output_prefix,
                                     LZW_FRACTION_CUTOFFS)
    # copy back to aws
    uploader_start(
        os.path.join(RESULT_DIR, LZW_OUT1), SAMPLE_S3_OUTPUT_PATH + "/")
    if len(input_fas) == 2:
        uploader_start(
            os.path.join(RESULT_DIR, LZW_OUT2), SAMPLE_S3_OUTPUT_PATH + "/")
    write_to_log("finished job")


def run_bowtie2(input_fas, uploader_start):
    # check if genome downloaded already
    genome_dir = fetch_genome(BOWTIE2_GENOME)

    # If we are here, we are also going to need a gsnap genome later;  start fetching it now.
    # This is actually THE PERFECT PLACE to initiate this fetch.  When we are running from
    # scratch, there is plenty of time to download the gsnap genome while bowtie is running.
    # When we are doing a lazy rerun, this function gets skipped, and the fetching of gsnap
    # genome is not initiated.  That's brilliant -- we don't fetch the gsnap genome if
    # we won't be needing it, and lazy reruns are very quick.
    threading.Thread(target=fetch_genome, args=[GSNAP_GENOME]).start()
    # the file structure looks like "bowtie2_genome/GRCh38.primary_assembly.genome.3.bt2"
    # the code below will handle up to "bowtie2_genome/GRCh38.primary_assembly.genome.99.bt2" but not 100
    local_genome_dir_ls = execute_command_with_output(
        "ls {genome_dir}/*.bt2*".format(genome_dir=genome_dir))
    genome_basename = local_genome_dir_ls.split("\n")[0][:-6]
    if genome_basename[-1] == '.':
        genome_basename = genome_basename[:-1]
    bowtie2_params = [
        BOWTIE2, '-q', '-p',
        str(multiprocessing.cpu_count()), '-x', genome_basename, '-f',
        '--very-sensitive-local', '-S', RESULT_DIR + '/' + BOWTIE2_OUT
    ]
    if len(input_fas) == 2:
        bowtie2_params.extend(['-1', input_fas[0], '-2', input_fas[1]])
    else:
        bowtie2_params.extend(['-U', input_fas[0]])
    execute_command(" ".join(bowtie2_params))
    write_to_log("finished alignment")
    # extract out unmapped files from sam
    output_prefix = RESULT_DIR + '/' + EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1[:
                                                                             -8]
    if len(input_fas) == 2:
        generate_unmapped_pairs_from_sam(RESULT_DIR + '/' + BOWTIE2_OUT,
                                         output_prefix)
    else:
        generate_unmapped_singles_from_sam(RESULT_DIR + '/' + BOWTIE2_OUT,
                                           output_prefix)
    uploader_start(
        os.path.join(RESULT_DIR, BOWTIE2_OUT), SAMPLE_S3_OUTPUT_PATH + "/")
    uploader_start(
        os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1),
        SAMPLE_S3_OUTPUT_PATH + "/")
    if len(input_fas) == 2:
        uploader_start(
            os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2),
            SAMPLE_S3_OUTPUT_PATH + "/")
        uploader_start(
            os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT3),
            SAMPLE_S3_OUTPUT_PATH + "/")
    write_to_log("extracted unmapped fragments from SAM file")


# Can remove this once the todo below, issue #173, is addressed.
class SkipGsnap(Exception):
    pass


def run_gsnap_filter(input_fas, uploader_start):
    # Unpack the gsnap genome
    genome_dir = fetch_genome(GSNAP_GENOME, strict=False)
    if genome_dir is None:
        # Apparently if the GSNAP_GENOME file doesn't exist, we are supposed to skip this step.
        # TODO (yunfang):  An independent way to specify whether this step should be executed,
        # so that operational errors don't just silently cause the step to be skipped. See #173.
        raise SkipGsnap()
    gsnap_base_dir = os.path.dirname(genome_dir)
    gsnap_index_name = os.path.basename(genome_dir)

    # Run Gsnap
    gsnap_params = [
        GSNAPL, '-A sam', '--batch=0', '--use-shared-memory=0',
        '--gmap-mode=all', '--npaths=1', '--ordered', '-t 32',
        '--max-mismatches=40', '-D', gsnap_base_dir, '-d', gsnap_index_name,
        '-o', RESULT_DIR + '/' + GSNAP_FILTER_SAM
    ]
    gsnap_params += input_fas
    execute_command(" ".join(gsnap_params))
    write_to_log("finished alignment")
    # extract out unmapped files from sam
    output_prefix = RESULT_DIR + '/' + EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1[:-8]
    if len(input_fas) == 2:
        generate_unmapped_pairs_from_sam(RESULT_DIR + '/' + GSNAP_FILTER_SAM,
                                         output_prefix)
    else:
        generate_unmapped_singles_from_sam(RESULT_DIR + '/' + GSNAP_FILTER_SAM,
                                           output_prefix)
    uploader_start(
        os.path.join(RESULT_DIR, GSNAP_FILTER_SAM),
        SAMPLE_S3_OUTPUT_PATH + "/")
    uploader_start(
        os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1),
        SAMPLE_S3_OUTPUT_PATH + "/")
    if len(input_fas) == 2:
        uploader_start(
            os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT2),
            SAMPLE_S3_OUTPUT_PATH + "/")
        uploader_start(
            os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT3),
            SAMPLE_S3_OUTPUT_PATH + "/")
    write_to_log("extracted unmapped fragments from SAM file for gsnap output")


@retry
def upload_with_retries(from_f, to_f):
    execute_command("aws s3 cp --quiet {from_f} {to_f}".format(
        from_f=from_f, to_f=to_f))


def upload(from_f, to_f, status, status_lock=threading.RLock()):
    try:
        with iostream_uploads:  # Limit concurrent uploads so as not to stall the pipeline.
            with iostream:  # Still counts toward the general semaphore.
                upload_with_retries(from_f, to_f)
            with status_lock:
                status[from_f] = "success"
    except:
        with status_lock:
            status[from_f] = "error"
        raise


def unzip_to_file(from_f, to_f):
    execute_command("gzip -dc {from_f} > {to_f}".format(
        from_f=from_f, to_f=to_f))


def run_host_filtering(fastq_files, initial_file_type_for_log, lazy_run, stats,
                       prefiltered):
    number_of_input_files = len(fastq_files)
    target_outputs = TARGET_OUTPUTS_SINGLE
    if number_of_input_files == 2:
        target_outputs = TARGET_OUTPUTS_PAIRED

    uploader_status = {}
    uploader_threads = []

    def uploader_start(from_f, to_f):
        t = threading.Thread(
            target=upload, args=[from_f, to_f, uploader_status])
        t.start()
        uploader_threads.append(t)

    def uploader_check_wait_all():
        for t in uploader_threads:
            t.join()
        for filename, status in uploader_status.iteritems():
            assert status == "success", "Bad upload status {} for file {}".format(
                status, filename)

    if prefiltered:
        # Move input in place of bowtie output (as it represents bowtie output from another run).
        btos = [
            EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1,
            EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2
        ]
        unzip_threads = []
        for i, fname in enumerate(fastq_files):
            assert fname.endswith(".fasta") or fname.endswith(
                ".fasta.gz"
            ), "Prefiltered input is not a fasta file: {fname}".format(
                fname=fname)
            if fname.endswith(".fasta"):
                execute_command("mv {fname} {bto}".format(
                    fname=fname, bto=os.path.join(RESULT_DIR, btos[i])))
            else:
                t = MyThread(
                    target=unzip_to_file,
                    args=[fname, os.path.join(RESULT_DIR, btos[i])])
                t.start()
                unzip_threads.append(t)
        for t in unzip_threads:
            t.join()
            assert t.completed and not t.exception
    else:
        # If not pre-filtered, run STAR
        # Getting total_reads and file truncation status from STAR
        total_counts_from_star = {}
        log_params = return_merged_dict(
            DEFAULT_LOG_PARAMS, {
                "title": "STAR",
                "version_file_s3": STAR_BOWTIE_VERSION_FILE_S3,
                "output_version_file": os.path.join(RESULT_DIR, VERSION_OUT)
            })
        run_and_log_s3(log_params, target_outputs["run_star"], lazy_run,
                       SAMPLE_S3_OUTPUT_PATH, run_star, fastq_files,
                       uploader_start, total_counts_from_star)
        if not total_counts_from_star.get('total_reads'):
            # Total reads not set. Most likely it's lazy run. Will have to actually count the reads.
            # TODO: Remove this when we also lazy load the stats.json file
            max_reads = MAX_INPUT_READS * len(fastq_files)
            total_reads = count_reads(fastq_files[0],
                                      initial_file_type_for_log, max_reads)
            if total_reads >= max_reads:
                total_reads = max_reads
                total_counts_from_star['truncated'] = 1
            total_counts_from_star['total_reads'] = total_reads

        stats.data.append(total_counts_from_star)
        stats.count_reads(
            "run_star",
            before_filename=fastq_files[0],
            before_filetype=initial_file_type_for_log,
            after_filename=os.path.join(RESULT_DIR, STAR_OUT1),
            after_filetype=initial_file_type_for_log)

        if total_counts_from_star.get('truncated'):
            # Upload the truncation file to notify web that the input files are truncated
            execute_command("echo %d | aws s3 cp - %s/%s" %
                            (total_counts_from_star['total_reads'],
                             SAMPLE_S3_OUTPUT_PATH, INPUT_TRUNCATED_FILE))

        # Run PriceSeqFilter
        log_params = return_merged_dict(DEFAULT_LOG_PARAMS,
                                        {"title": "PriceSeqFilter"})
        input_files = [os.path.join(RESULT_DIR, STAR_OUT1)]
        if number_of_input_files == 2:
            input_files.append(os.path.join(RESULT_DIR, STAR_OUT2))

        run_and_log_s3(log_params, target_outputs["run_priceseqfilter"],
                       lazy_run, SAMPLE_S3_OUTPUT_PATH, run_priceseqfilter,
                       input_files, uploader_start)
        stats.count_reads(
            "run_priceseqfilter",
            before_filename=os.path.join(RESULT_DIR, STAR_OUT1),
            before_filetype=initial_file_type_for_log,
            after_filename=os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1),
            after_filetype=initial_file_type_for_log)

        # Run FASTQ to FASTA
        if "fastq" in FILE_TYPE:
            log_params = return_merged_dict(DEFAULT_LOG_PARAMS,
                                            {"title": "FASTQ to FASTA"})
            input_files = [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1)]
            next_inputs = [os.path.join(RESULT_DIR, FQ2FA_OUT1)]
            if number_of_input_files == 2:
                input_files.append(
                    os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2))
                next_inputs.append(os.path.join(RESULT_DIR, FQ2FA_OUT2))
            run_and_log_s3(log_params, target_outputs["run_fq2fa"], lazy_run,
                           SAMPLE_S3_OUTPUT_PATH, run_fq2fa, input_files,
                           uploader_start)
        else:
            next_inputs = [os.path.join(RESULT_DIR, PRICESEQFILTER_OUT1)]
            if number_of_input_files == 2:
                next_inputs.append(
                    os.path.join(RESULT_DIR, PRICESEQFILTER_OUT2))

        # Run CD-HIT-DUP
        log_params = return_merged_dict(DEFAULT_LOG_PARAMS,
                                        {"title": "CD-HIT-DUP"})
        run_and_log_s3(log_params, target_outputs["run_cdhitdup"], lazy_run,
                       SAMPLE_S3_OUTPUT_PATH, run_cdhitdup, next_inputs,
                       uploader_start)
        stats.count_reads(
            "run_cdhitdup",
            before_filename=next_inputs[0],
            before_filetype="fasta_paired",
            after_filename=os.path.join(RESULT_DIR, CDHITDUP_OUT1),
            after_filetype="fasta_paired")

        # Run LZW filter
        log_params = return_merged_dict(DEFAULT_LOG_PARAMS,
                                        {"title": "LZW filter"})
        input_files = [os.path.join(RESULT_DIR, CDHITDUP_OUT1)]
        if number_of_input_files == 2:
            input_files.append(os.path.join(RESULT_DIR, CDHITDUP_OUT2))
        run_and_log_s3(log_params, target_outputs["run_lzw"], lazy_run,
                       SAMPLE_S3_OUTPUT_PATH, run_lzw, input_files,
                       uploader_start)
        stats.count_reads(
            "run_lzw",
            before_filename=os.path.join(RESULT_DIR, CDHITDUP_OUT1),
            before_filetype="fasta_paired",
            after_filename=os.path.join(RESULT_DIR, LZW_OUT1),
            after_filetype="fasta_paired")

        # Run Bowtie
        log_params = return_merged_dict(DEFAULT_LOG_PARAMS,
                                        {"title": "bowtie2"})
        input_files = [os.path.join(RESULT_DIR, LZW_OUT1)]
        if number_of_input_files == 2:
            input_files.append(os.path.join(RESULT_DIR, LZW_OUT2))

        # Upload Bowtie2 results
        run_and_log_s3(log_params, target_outputs["run_bowtie2"], lazy_run,
                       SAMPLE_S3_OUTPUT_PATH, run_bowtie2, input_files,
                       uploader_start)
        stats.count_reads(
            "run_bowtie2",
            before_filename=os.path.join(RESULT_DIR, LZW_OUT1),
            before_filetype="fasta_paired",
            after_filename=os.path.join(RESULT_DIR,
                                        EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1),
            after_filetype="fasta_paired")

    # Run GSNAP against host genomes (only available for Human as of 5/1/2018)
    # GSNAP may run again even for pre-filtered inputs
    try:
        input_files = [
            os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1)
        ]
        if number_of_input_files == 2:
            input_files.append(
                os.path.join(RESULT_DIR,
                             EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2))

        # Skip GSNAP if the number of reads is too big
        # TODO: move gsnap filter to after subsampling
        if stats.data[-1]['reads_after'] > MAX_GNSAP_FILTER_READS:
            raise SkipGsnap()
        log_params = return_merged_dict(DEFAULT_LOG_PARAMS,
                                        {"title": "run_gsnap_filter"})
        run_and_log_s3(log_params, target_outputs["run_gsnap_filter"],
                       lazy_run, SAMPLE_S3_OUTPUT_PATH, run_gsnap_filter,
                       input_files, uploader_start)
        stats.count_reads(
            "run_gsnap_filter",
            before_filename=os.path.join(
                RESULT_DIR, EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1),
            before_filetype="fasta_paired",
            after_filename=os.path.join(RESULT_DIR,
                                        EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1),
            after_filetype="fasta_paired")
    except SkipGsnap:
        write_to_log(
            "Skipping gsnap is for prefilterd input or too many reads")
        pass

    # Finalize the remaining reads
    stats.set_remaining_reads()
    uploader_check_wait_all()


def upload_pipeline_version_file():
    execute_command("echo %s > %s" % (PIPELINE_VERSION, PIPELINE_VERSION_OUT))
    execute_command("aws s3 cp %s %s/" % (PIPELINE_VERSION_OUT, OUTPUT_BUCKET))


def run_stage1(lazy_run=True):
    execute_command("mkdir -p %s %s %s %s" % (SAMPLE_DIR, FASTQ_DIR,
                                              RESULT_DIR, SCRATCH_DIR))
    execute_command("mkdir -p %s " % REF_DIR)

    # configure logger
    log_file = "%s/%s.%s.txt" % (RESULT_DIR, LOGS_OUT_BASENAME,
                                 AWS_BATCH_JOB_ID)
    configure_logger(log_file)
    print("Starting stage...")

    # Get list of input files to fetch
    command = "aws s3 ls %s/ | grep '\\.%s$'" % (SAMPLE_S3_INPUT_PATH,
                                                 FILE_TYPE)
    output = execute_command_with_output(command).rstrip().split("\n")
    input_fetch_threads = []

    def fetch_input(input_basename):
        fetch_from_s3(
            os.path.join(SAMPLE_S3_INPUT_PATH, input_basename),
            FASTQ_DIR,
            allow_s3mi=True,
            auto_unzip=False)

    # Fetch input files with multiple threads
    for line in output:
        m = re.match(".*?([^ ]*." + re.escape(FILE_TYPE) + ")", line)
        if m:
            t = MyThread(target=fetch_input, args=[m.group(1)])
            t.start()
            input_fetch_threads.append(t)
        else:
            print("%s doesn't match %s" % (line, FILE_TYPE))
    for t in input_fetch_threads:
        # Check thread completion
        t.join()
        assert t.completed and not t.exception

    # Check FASTQ files
    fastq_files = execute_command_with_output(
        "ls %s/*.%s" % (FASTQ_DIR, FILE_TYPE)).rstrip().split("\n")
    if len(fastq_files) not in [1, 2]:
        write_to_log(
            "Number of input files was neither 1 nor 2. Aborting computation.")
        return  # only support either 1 file or 2 (paired) files

    initial_file_type_for_log = "fasta"
    if "fastq" in FILE_TYPE:
        initial_file_type_for_log = "fastq"
    if len(fastq_files) == 2:
        initial_file_type_for_log += "_paired"

    # Instantiate a stats instance
    stats = StatsFile(STATS_OUT, RESULT_DIR, None, SAMPLE_S3_OUTPUT_PATH)

    # Download total_reads.json input, if present.  This is only provided with post-filtered inputs,
    # where we don't have the reads prior to host filtering.
    # Record total number of input reads
    try:
        stats_in = StatsFile(STATS_IN, RESULT_DIR, SAMPLE_S3_INPUT_PATH,
                             SAMPLE_S3_OUTPUT_PATH)
        stats_in.load_from_s3()
        total_reads = stats_in.get_total_reads()
        assert total_reads == int(total_reads)
        write_to_log(
            "Post-filtered input with {total_reads} original total reads.".
            format(total_reads=total_reads))
    except:
        total_reads = None
        stats_in = None
        write_to_log("Unfiltered input. Need host filtering")

    if total_reads is not None:  # set total reads if available
        stats.data.append({'total_reads': total_reads})

    # Run host filtering
    run_host_filtering(fastq_files, initial_file_type_for_log, lazy_run, stats,
                       stats_in is not None)

    stats.save_to_s3()

    write_to_log("Host filtering complete")
    upload_log_file(SAMPLE_S3_OUTPUT_PATH)
