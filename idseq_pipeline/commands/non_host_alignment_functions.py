import os
import subprocess
import json
import shelve
import time
import math
import threading
import shutil
import traceback
import random
from collections import defaultdict
from .common import *  #pylint: disable=wildcard-import

# Settings per job. By default jobs are subsampled to <= 100 chunks.
MAX_CHUNKS_IN_FLIGHT_GSNAP = 32
MAX_CHUNKS_IN_FLIGHT_RAPSEARCH = 32
chunks_in_flight_gsnap = threading.Semaphore(MAX_CHUNKS_IN_FLIGHT_GSNAP)
chunks_in_flight_rapsearch = threading.Semaphore(
    MAX_CHUNKS_IN_FLIGHT_RAPSEARCH)

# Dispatch at most this many chunks per minute, to ensure fairness
# amongst jobs regardless of job size (not the best way to do it,
# just for now).
MAX_DISPATCHES_PER_MINUTE = 10

# poll this many random servers per wait_for_instance_ip
MAX_INSTANCES_TO_POLL = 8

MAX_POLLING_LATENCY = 10  # seconds

# If no instance is available, we should refresh our list of instances,
# to pick up new instances added by auto-scaling. Wait at least this long
# between refreshes to stay within AWS account limits.
MIN_INTERVAL_BETWEEN_DESCRIBE_INSTANCES = 180

# Refresh at least every 30 minutes
MAX_INTERVAL_BETWEEN_DESCRIBE_INSTANCES = 900

# data directories
ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + '/idseq/data'  # generated data go here
REF_DIR = ROOT_DIR + '/idseq/ref'  # referene genome / ref databases go here
# input files
EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta'
EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta'
EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT3 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.merged.fasta'

EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1 = 'unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta'
EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT2 = 'unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta'
EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT3 = 'unmapped.gsnap_filter.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.merged.fasta'
# output files
UNIDENTIFIED_FASTA_OUT = 'unidentified.fasta'
DEPRECATED_BOOBYTRAPPED_COMBINED_JSON_OUT = 'idseq_web_sample.json'
LOGS_OUT_BASENAME = 'log'
STATS_OUT = 'stats.json'
VERSION_OUT = 'versions.json'
MULTIHIT_GSNAPL_OUT = 'multihit.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
SUMMARY_MULTIHIT_GSNAPL_OUT = 'summary.multihit.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.tab'
DEDUP_MULTIHIT_GSNAPL_OUT = 'dedup.multihit.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
MULTIHIT_NT_JSON_OUT = 'nt_multihit_idseq_web_sample.json'
MULTIHIT_RAPSEARCH_OUT = 'multihit.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
SUMMARY_MULTIHIT_RAPSEARCH_OUT = 'summary.multihit.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.tab'
DEDUP_MULTIHIT_RAPSEARCH_OUT = 'dedup.multihit.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
MULTIHIT_NR_JSON_OUT = 'nr_multihit_idseq_web_sample.json'
MULTIHIT_COMBINED_JSON_OUT = 'multihit_idseq_web_sample.json'
ACCESSION_ANNOTATED_FASTA = 'accessions.rapsearch2.gsnapl.fasta'

# arguments from environment variables
SKIP_DEUTERO_FILTER = int(os.environ.get('SKIP_DEUTERO_FILTER', 0))
SUBSAMPLE = os.environ.get(
    'SUBSAMPLE'
)  # number of read pairs to subsample to, before gsnap/rapsearch
FASTQ_BUCKET = os.environ.get('FASTQ_BUCKET')
INPUT_BUCKET = os.environ.get('INPUT_BUCKET')
FILE_TYPE = os.environ.get('FILE_TYPE', 'fastq.gz')
OUTPUT_BUCKET = os.environ.get('OUTPUT_BUCKET')
AWS_BATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID', 'local')
ENVIRONMENT = os.environ.get('ENVIRONMENT', 'prod')
SAMPLE_S3_FASTQ_PATH = FASTQ_BUCKET.rstrip('/')
SAMPLE_S3_INPUT_PATH = INPUT_BUCKET.rstrip('/')
SAMPLE_S3_OUTPUT_PATH = OUTPUT_BUCKET.rstrip('/')
SAMPLE_S3_OUTPUT_CHUNKS_PATH = os.path.join(SAMPLE_S3_OUTPUT_PATH, "chunks")
SAMPLE_NAME = SAMPLE_S3_INPUT_PATH[5:].rstrip('/').replace('/', '-')
SAMPLE_DIR = DEST_DIR + '/' + SAMPLE_NAME
FASTQ_DIR = SAMPLE_DIR + '/fastqs'
RESULT_DIR = SAMPLE_DIR + '/results'


def result_dir(basename):
    return os.path.join(RESULT_DIR, basename)


CHUNKS_RESULT_DIR = result_dir("chunks")
DEFAULT_LOG_PARAMS = {"sample_s3_output_path": SAMPLE_S3_OUTPUT_PATH}

# For reproducibility of random operations
random.seed(hash(SAMPLE_NAME))

# versioning
## For now, index updates are infrequent and we can get their versions from S3.
## If updates ever become frequent, we may want to check instead which version is actually on the
## machine taking the job, possibly out of sync with the newest version in S3.

base_s3 = 's3://idseq-database/alignment_indexes'
base_dt = '2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime'
GSNAP_VERSION_FILE_S3 = ("%s/%s/nt_k16.version.txt" % (base_s3, base_dt))
RAPSEARCH_VERSION_FILE_S3 = ("%s/%s/nr_rapsearch.version.txt" % (base_s3,
                                                                 base_dt))

# compute capacity
GSNAPL_MAX_CONCURRENT = 3  # number of gsnapl jobs allowed to run concurrently on 1 machine
RAPSEARCH2_MAX_CONCURRENT = 6
GSNAPL_CHUNK_SIZE = 15000  # number of fasta records in a chunk so it runs in ~10 minutes on i3.16xlarge
RAPSEARCH_CHUNK_SIZE = 10000

# references
# from common import ACCESSION2TAXID
base_s3 = 's3://idseq-database/taxonomy'
base_dt = '2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime'
DEUTEROSTOME_TAXIDS = ("%s/%s/deuterostome_taxids.txt" % (base_s3, base_dt))

# from common import LINEAGE_SHELF


# convenience functions
def count_lines(input_file):
    return int(
        execute_command_with_output(
            "wc -l %s" % input_file).strip().split()[0])


def count_lines_in_paired_files(input_files):
    distinct_counts = list(set(map(count_lines, input_files)))
    assert len(
        distinct_counts
    ) == 1, "Mismatched line counts in supposedly paired files: {}".format(
        input_files)
    return distinct_counts[0]


def subsample_helper(input_file, records_to_keep, type_, output_file):
    """Subsample the FASTA file by either read IDs or record indices to save.
    record_indices are used first to sample from all the records in the file.
    The paired FASTA case uses read_ids to subsample the merged file.
    """

    record_number = 0
    msg = "records_to_keep is not a set."
    assert isinstance(records_to_keep, set), msg
    kept_read_ids = set()
    kept_reads_count = 0
    sequence_basename = ""

    with open(input_file, 'rb') as input_f:
        with open(output_file, 'wb') as output_f:
            # Iterate through the FASTA file records
            sequence_name = input_f.readline()
            sequence_data = input_f.readline()

            while len(sequence_name) > 0 and len(sequence_data) > 0:
                if type_ == "read_ids":
                    sequence_basename = sequence_name.rstrip().rsplit('/',
                                                                      1)[0]
                    condition = sequence_basename in records_to_keep
                elif type_ == "record_indices":
                    condition = record_number in records_to_keep
                    sequence_basename = sequence_name.rstrip()
                else:
                    condition = True

                if condition:  # We should keep this entry
                    output_f.write(sequence_name)
                    output_f.write(sequence_data)
                    kept_reads_count += 1
                    kept_read_ids.add(sequence_basename)
                sequence_name = input_f.readline()
                sequence_data = input_f.readline()
                record_number += 1

    # Error checking for mismatch
    if type_ == "record_indices" and len(kept_read_ids) != len(records_to_keep):
        msg = "WARNING: Repeated read IDs in {input_file}: Found {kri} read " \
              "IDs in {rtk} sampled rows.".format(input_file=input_file,
                                                  kri=len(kept_read_ids),
                                                  rtk=len(records_to_keep))
        write_to_log(msg)

    if type_ == "read_ids" and kept_read_ids != records_to_keep:
        missing = records_to_keep - kept_read_ids
        examples = sorted(random.sample(missing, min(10, len(missing))))
        msg = "Not all desired read IDs were found in the file: {}\nMissing: {}".format(
            input_file, examples)
        assert kept_read_ids == records_to_keep, msg
    return kept_read_ids, kept_reads_count


# Note dedicated random stream for just this function, so that arbitrary
# other use of randomness (e.g. for I/O retry delays) will not perturb the
# sub-sampling stream
def subsample_single_fasta(input_files_basename,
                           target_n_reads,
                           randgen=random.Random(x=hash(SAMPLE_NAME))):
    """Randomly subsample a single unpaired FASTA file."""
    input_file = result_dir(input_files_basename)
    # Each FASTA record spans 2 lines
    total_records = count_lines(input_file) // 2
    write_to_log("total unpaired reads: %d" % total_records)
    write_to_log("target reads: %d" % target_n_reads)
    if total_records <= target_n_reads:
        return input_file

    # Randomly select a subset to process
    records_to_keep = set(
        randgen.sample(xrange(total_records), target_n_reads))
    subsample_prefix = "subsample_%d" % target_n_reads

    input_dir, input_basename = os.path.split(input_file)
    output_basename = "%s.%s" % (subsample_prefix, input_basename)
    output_file = os.path.join(input_dir, output_basename)
    subsample_helper(input_file, records_to_keep, "record_indices",
                     output_file)
    return output_file


def subsample_paired_fastas(input_files_basenames,
                            merged_file_basename,
                            target_n_reads,
                            randgen=random.Random(x=hash(SAMPLE_NAME))):
    """Randomly subsample a paired FASTA file. A paired read means reading
    from both ends of a DNA/RNA fragment.
    """
    input_files = [result_dir(f) for f in input_files_basenames]
    merged_file = result_dir(merged_file_basename)
    # Each FASTA record spans 2 lines
    total_records = count_lines_in_paired_files(input_files) // 2
    write_to_log("total read pairs: %d" % total_records)
    write_to_log("target read pairs: %d" % target_n_reads)
    # Note: target_n_reads and total_records here really refer to numbers of
    # read PAIRS.
    if total_records <= target_n_reads:
        return input_files_basenames, merged_file_basename

    subsample_prefix = "subsample_%d" % target_n_reads
    records_to_keep = set(
        randgen.sample(xrange(total_records), target_n_reads))
    subsampled_files = []
    known_kept_read_ids = set()

    # Subsample the paired files and record read IDs kept
    for input_file in input_files:
        input_dir = os.path.split(input_file)[0]
        input_basename = os.path.split(input_file)[1]
        output_basename = "%s.%s" % (subsample_prefix, input_basename)
        output_file = os.path.join(input_dir, output_basename)

        kept_read_ids, kept_reads_count = subsample_helper(
            input_file, records_to_keep, "record_indices", output_file)
        assert kept_reads_count == target_n_reads
        subsampled_files.append(output_basename)
        known_kept_read_ids |= kept_read_ids  # Set union

    # Subsample the merged file to the same read IDs
    input_dir = os.path.split(merged_file)[0]
    input_basename = os.path.split(merged_file)[1]
    subsampled_merged_file = "%s.%s" % (subsample_prefix, input_basename)
    output_file = os.path.join(input_dir, subsampled_merged_file)
    _, kept_reads_count = subsample_helper(merged_file, known_kept_read_ids,
                                           "read_ids", output_file)
    num_desired = target_n_reads * len(input_files_basenames)

    if kept_reads_count != num_desired:
        # TODO: Surface these warnings in the webapp.
        msg = "ERROR: Improperly paired reads. Total reads in sampled merged " \
              "fasta: {krc}".format(krc=kept_reads_count)
        write_to_log(msg)
    assert kept_reads_count == num_desired
    return subsampled_files, subsampled_merged_file


def concatenate_files(file_list, output_file):
    with open(output_file, 'wb') as outf:
        for f in file_list:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, outf)


def log_corrupt(m8_file, line):
    msg = m8_file + " is corrupt at line:\n" + line + "\n----> delete it and its corrupt ancestors before restarting run"
    write_to_log(msg)
    return msg


def iterate_m8(m8_file, debug_caller=None, logging_interval=25000000):
    """Generate an iterator over the m8 file and return (read_id,
    accession_id, percent_id, alignment_length, e_value, line) for each line.
    Work around and warn about any invalid hits detected.
    """
    invalid_hits = 0
    last_invalid_line = None
    with open(m8_file, 'rb') as m8f:
        line_count = 0
        for line in m8f:
            line_count += 1
            if line and line[0] == '#':
                # skip comment lines
                continue
            parts = line.split("\t")
            # Must have at least 12 parts per line
            assert len(parts) >= 12, log_corrupt(m8_file, line)

            read_id = parts[0]
            accession_id = parts[1]
            percent_id = float(parts[2])
            alignment_length = int(parts[3])
            e_value = float(parts[10])

            # GSNAP outputs bogus alignments (non-positive length /
            # impossible percent identity / NaN e-value) sometimes,
            # and usually they are not the only assignment, so rather than
            # killing the job, we just skip them. If we don't filter these
            # out here, they will override the good data when computing min(
            # evalue), pollute averages computed in the json, and cause the
            # webapp loader to crash as the Rails JSON parser cannot handle
            # NaNs. Test if e_value != e_value to test if e_value is NaN
            # because NaN != NaN.
            if alignment_length <= 0 or not -0.25 < percent_id < 100.25 or e_value != e_value:
                invalid_hits += 1
                last_invalid_line = line
                continue
            if debug_caller and line_count % logging_interval == 0:
                msg = "Scanned {} m8 lines from {} for {}, and going.".format(
                    line_count, m8_file, debug_caller)
                write_to_log(msg)
            yield (read_id, accession_id, percent_id, alignment_length,
                   e_value, line)

    # Warn about any invalid hits outputted by GSNAP
    if invalid_hits:
        msg = "Found {} invalid hits in {};  last invalid hit line: {}".format(
            invalid_hits, m8_file, last_invalid_line)
        write_to_log(msg, warning=True)
    if debug_caller:
        msg = "Scanned all {} m8 lines from {} for {}.".format(
            line_count, m8_file, debug_caller)
        write_to_log(msg)


def annotate_fasta_with_accessions(merged_input_fasta, nt_m8, nr_m8,
                                   output_fasta):
    def get_map(m8_file):
        return dict((read_id, accession_id)
                    for read_id, accession_id, _percent_id,
                    _alignment_length, _e_value, _line in iterate_m8(
                        m8_file, "annotate_fasta_with_accessions"))

    nt_map = get_map(nt_m8)
    nr_map = get_map(nr_m8)

    with open(merged_input_fasta, 'rb') as input_fasta_f:
        with open(output_fasta, 'wb') as output_fasta_f:
            sequence_name = input_fasta_f.readline()
            sequence_data = input_fasta_f.readline()

            while sequence_name and sequence_data:
                read_id = sequence_name.rstrip().lstrip('>')
                # Need to annotate NR then NT in this order for alignment viz
                new_read_name = "NR:{nr_accession}:NT:{nt_accession}:{read_id}".format(
                    nr_accession=nr_map.get(read_id, ''),
                    nt_accession=nt_map.get(read_id, ''),
                    read_id=read_id)
                output_fasta_f.write(">%s\n" % new_read_name)
                output_fasta_f.write(sequence_data)
                sequence_name = input_fasta_f.readline()
                sequence_data = input_fasta_f.readline()

    execute_command("aws s3 cp --quiet %s %s/" % (output_fasta,
                                                  SAMPLE_S3_OUTPUT_PATH))


@run_in_subprocess
def generate_taxon_count_json_from_m8(
        m8_file, hit_level_file, e_value_type, count_type, lineage_map_path,
        deuterostome_path, total_reads, remaining_reads, output_file):
    # Parse through hit file and m8 input file and format a JSON file with
    # our desired attributes, including aggregated statistics.

    if not SKIP_DEUTERO_FILTER:
        taxids_to_remove = read_file_into_set(deuterostome_path)

    def any_hits_to_remove(hits):
        if SKIP_DEUTERO_FILTER:
            return False
        for taxid in hits:
            if int(taxid) >= 0 and taxid in taxids_to_remove:
                return True
        return False

    # Setup
    aggregation = {}
    hit_f = open(hit_level_file, 'rb')
    m8_f = open(m8_file, 'rb')
    # Lines in m8_file and hit_level_file correspond (same read_id)
    hit_line = hit_f.readline()
    m8_line = m8_f.readline()
    lineage_map = shelve.open(lineage_map_path)
    NUM_RANKS = len(NULL_LINEAGE)
    # See https://en.wikipedia.org/wiki/Double-precision_floating-point_format
    MIN_NORMAL_POSITIVE_DOUBLE = 2.0**-1022

    while hit_line and m8_line:
        # Retrieve data values from files
        hit_line_columns = hit_line.rstrip("\n").split("\t")
        _read_id = hit_line_columns[0]
        hit_level = hit_line_columns[1]
        hit_taxid = hit_line_columns[2]
        if int(hit_level) < 0:  # Skip negative levels and continue
            hit_line = hit_f.readline()
            m8_line = m8_f.readline()
            continue

        # m8 files correspond to BLAST tabular output format 6:
        # Columns: read_id | _ref_id | percent_identity | alignment_length...
        #
        # * read_id = query (e.g., gene) sequence id
        # * _ref_id = subject (e.g., reference genome) sequence id
        # * percent_identity = percentage of identical matches
        # * alignment_length = length of the alignments
        # * e_value = the expect value
        #
        # See:
        # * http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
        # * http://www.metagenomics.wiki/tools/blast/evalue

        m8_line_columns = m8_line.split("\t")
        msg = "read_ids in %s and %s do not match: %s vs. %s" % (
            os.path.basename(m8_file), os.path.basename(hit_level_file),
            m8_line_columns[0], hit_line_columns[0])
        assert m8_line_columns[0] == hit_line_columns[0], msg
        percent_identity = float(m8_line_columns[2])
        alignment_length = float(m8_line_columns[3])
        e_value = float(m8_line_columns[10])

        # These have been filtered out before the creation of m8_f and hit_f
        assert alignment_length > 0
        assert -0.25 < percent_identity < 100.25
        assert e_value == e_value
        if e_value_type != 'log10':
            # Subtle, but this excludes positive subnormal numbers.
            assert MIN_NORMAL_POSITIVE_DOUBLE <= e_value
            e_value = math.log10(e_value)

        # Retrieve the taxon lineage and mark meaningless calls with fake
        # taxids.
        hit_taxids_all_levels = lineage_map.get(hit_taxid, NULL_LINEAGE)
        cleaned_hit_taxids_all_levels = validate_taxid_lineage(
            hit_taxids_all_levels, hit_taxid, hit_level)
        assert NUM_RANKS == len(cleaned_hit_taxids_all_levels)

        if not any_hits_to_remove(cleaned_hit_taxids_all_levels):
            # Aggregate each level and collect statistics
            agg_key = tuple(cleaned_hit_taxids_all_levels)
            while agg_key:
                agg_bucket = aggregation.get(agg_key)
                if not agg_bucket:
                    agg_bucket = {
                        'count': 0,
                        'sum_percent_identity': 0.0,
                        'sum_alignment_length': 0.0,
                        'sum_e_value': 0.0
                    }
                    aggregation[agg_key] = agg_bucket
                agg_bucket['count'] += 1
                agg_bucket['sum_percent_identity'] += percent_identity
                agg_bucket['sum_alignment_length'] += alignment_length
                agg_bucket['sum_e_value'] += e_value
                # Chomp off the lowest rank as we aggregate up the tree
                agg_key = agg_key[1:]

        hit_line = hit_f.readline()
        m8_line = m8_f.readline()

    # Produce the final output
    taxon_counts_attributes = []
    for agg_key, agg_bucket in aggregation.iteritems():
        count = agg_bucket['count']
        tax_level = NUM_RANKS - len(agg_key) + 1
        # TODO: Extend taxonomic ranks as indicated on the commented out lines.
        taxon_counts_attributes.append({
            "tax_id":
            agg_key[0],
            "tax_level":
            tax_level,
            # 'species_taxid' : agg_key[tax_level - 1] if tax_level == 1 else "-100",
            'genus_taxid':
            agg_key[2 - tax_level] if tax_level <= 2 else "-200",
            'family_taxid':
            agg_key[3 - tax_level] if tax_level <= 3 else "-300",
            # 'order_taxid' : agg_key[4 - tax_level] if tax_level <= 4 else "-400",
            # 'class_taxid' : agg_key[5 - tax_level] if tax_level <= 5 else "-500",
            # 'phyllum_taxid' : agg_key[6 - tax_level] if tax_level <= 6 else "-600",
            # 'kingdom_taxid' : agg_key[7 - tax_level] if tax_level <= 7 else "-700",
            # 'domain_taxid' : agg_key[8 - tax_level] if tax_level <= 8 else "-800",
            "count":
            count,
            "percent_identity":
            agg_bucket['sum_percent_identity'] / count,
            "alignment_length":
            agg_bucket['sum_alignment_length'] / count,
            "e_value":
            agg_bucket['sum_e_value'] / count,
            "count_type":
            count_type
        })
    output_dict = {
        "pipeline_output": {
            "total_reads": total_reads,
            "remaining_reads": remaining_reads,
            "taxon_counts_attributes": taxon_counts_attributes
        }
    }

    with open(output_file, 'wb') as outf:
        json.dump(output_dict, outf)
    execute_command("aws s3 cp --quiet %s %s/" % (output_file,
                                                  SAMPLE_S3_OUTPUT_PATH))


def combine_pipeline_output_json(inputPath1, inputPath2, outputPath, stats):
    total_reads = stats.get_total_reads()
    remaining_reads = stats.get_remaining_reads()
    with open(inputPath1) as inf1:
        input1 = json.load(inf1).get("pipeline_output")
    with open(inputPath2) as inf2:
        input2 = json.load(inf2).get("pipeline_output")
    taxon_counts_attributes = (input1.get("taxon_counts_attributes") +
                               input2.get("taxon_counts_attributes"))
    pipeline_output_dict = {
        "total_reads": total_reads,
        "remaining_reads": remaining_reads,
        "taxon_counts_attributes": taxon_counts_attributes
    }
    output_dict = {"pipeline_output": pipeline_output_dict}
    with open(outputPath, 'wb') as outf:
        json.dump(output_dict, outf)
    execute_command("aws s3 cp --quiet %s %s/" % (outputPath,
                                                  SAMPLE_S3_OUTPUT_PATH))


def read_file_into_set(file_name):
    with open(file_name) as f:
        S = set(x.rstrip() for x in f)
    S.discard('')
    return S


def environment_for_aligners(_environment):
    if _environment in ['development', 'staging', 'prod']:
        return 'prod'
    return ''


def fetch_key(environment, mutex=threading.RLock()):
    with mutex:
        key_s3_path = "s3://idseq-secrets/idseq-%s.pem" % environment_for_aligners(
            environment)
        key_name = os.path.basename(key_s3_path)
        key_path = REF_DIR + '/' + key_name
        if not os.path.exists(key_path):
            execute_command("aws s3 cp --quiet %s %s/" % (key_s3_path,
                                                          REF_DIR))
            execute_command("chmod 400 %s" % key_path)
        return key_path


@retry
def get_server_ips_work(service_name, environment):
    tag = "service"
    value = "%s-%s" % (service_name, environment_for_aligners(environment))
    describe_json = json.loads(
        execute_command_with_output(
            "aws ec2 describe-instances --filters 'Name=tag:%s,Values=%s' 'Name=instance-state-name,Values=running'"
            % (tag, value)))
    server_ips = []
    for reservation in describe_json["Reservations"]:
        for instance in reservation["Instances"]:
            server_ips += [
                instance["NetworkInterfaces"][0]["Association"]["PublicIp"]
            ]
    return server_ips


def get_server_ips(service_name,
                   environment,
                   aggressive=False,
                   cache={},
                   mutex=threading.RLock()):  #pylint: disable=dangerous-default-value
    try:
        with mutex:
            if aggressive:
                period = MIN_INTERVAL_BETWEEN_DESCRIBE_INSTANCES
            else:
                period = MAX_INTERVAL_BETWEEN_DESCRIBE_INSTANCES
            now = time.time()
            cache_key = (service_name, environment)
            if cache_key not in cache or now - cache[cache_key][0] >= period:
                # this may raise an exception when the AWS account rate limit is exceeded due to many concurrent jobs
                cache[cache_key] = (now,
                                    get_server_ips_work(
                                        service_name, environment))
            return cache[cache_key][1]
    except:
        # return [] causes a sleep of wait_seconds before retrying (see below)
        traceback.print_exc()
        return []


def wait_for_server_ip_work(service_name,
                            key_path,
                            remote_username,
                            environment,
                            max_concurrent,
                            chunk_id,
                            had_to_wait=[False]):  #pylint: disable=dangerous-default-value
    while True:
        write_to_log(
            "Chunk {chunk_id} of {service_name} is at third gate".format(
                chunk_id=chunk_id, service_name=service_name))
        instance_ips = get_server_ips(
            service_name, environment, aggressive=had_to_wait[0])
        instance_ips = random.sample(instance_ips,
                                     min(MAX_INSTANCES_TO_POLL,
                                         len(instance_ips)))
        ip_nproc_dict = {}
        dict_mutex = threading.RLock()
        dict_writable = True

        def poll_server(ip):
            # ServerAliveInterval to fix issue with containers keeping open
            # an SSH connection even after worker machines had finished
            # running.
            command = 'ssh -o "StrictHostKeyChecking no" -o "ConnectTimeout 5" -o "ServerAliveInterval 60" -i %s %s@%s "ps aux | grep %s | grep -v bash" || echo "error"' % (
                key_path, remote_username, ip, service_name)
            output = execute_command_with_output(
                command, timeout=MAX_POLLING_LATENCY).rstrip().split("\n")
            if output != ["error"]:
                with dict_mutex:
                    if dict_writable:
                        ip_nproc_dict[ip] = len(output) - 1

        poller_threads = []
        for ip in instance_ips:
            pt = threading.Thread(target=poll_server, args=[ip])
            pt.start()
            poller_threads.append(pt)
        for pt in poller_threads:
            pt.join(MAX_POLLING_LATENCY)
        with dict_mutex:
            # Any zombie threads won't be allowed to modify the dict.
            dict_writable = False
        if len(ip_nproc_dict) < len(poller_threads):
            write_to_log(
                "Only {} out of {} instances responded to polling;  {} threads are timing out.".
                format(
                    len(ip_nproc_dict), len(poller_threads),
                    len([pt for pt in poller_threads if pt.isAlive()])))
        write_to_log(
            "Chunk {chunk_id} of {service_name} is at fourth gate".format(
                chunk_id=chunk_id, service_name=service_name))
        if not ip_nproc_dict:
            have_capacity = False
        else:
            min_nproc_ip = min(ip_nproc_dict, key=ip_nproc_dict.get)
            min_nproc = ip_nproc_dict[min_nproc_ip]
            have_capacity = (min_nproc < max_concurrent)
        if have_capacity:
            had_to_wait[0] = False
            # Make an urn where each ip occurs more times if it has more free slots, and not at all if it lacks free slots
            urn = []
            for ip, nproc in ip_nproc_dict.iteritems():
                free_slots = max_concurrent - nproc
                if free_slots > 0:
                    weight = 2**free_slots - 1
                    urn.extend([ip] * weight)
            min_nproc_ip = random.choice(urn)
            free_slots = max_concurrent - ip_nproc_dict[min_nproc_ip]
            write_to_log("%s server %s has capacity %d. Kicking off " %
                         (service_name, min_nproc_ip, free_slots))
            return min_nproc_ip
        else:
            had_to_wait[0] = True
            wait_seconds = random.randint(
                max(20, MAX_POLLING_LATENCY), max(60, MAX_POLLING_LATENCY))
            write_to_log("%s servers busy. Wait for %d seconds" %
                         (service_name, wait_seconds))
            time.sleep(wait_seconds)


def wait_for_server_ip(service_name,
                       key_path,
                       remote_username,
                       environment,
                       max_concurrent,
                       chunk_id,
                       mutex=threading.RLock(),
                       mutexes={},
                       last_checks={}):  #pylint: disable=dangerous-default-value
    # We rate limit these to ensure fairness across jobs regardless of job size
    with mutex:
        if service_name not in mutexes:
            mutexes[service_name] = threading.RLock()
            last_checks[service_name] = [None]
        lc = last_checks[service_name]
        mx = mutexes[service_name]
    write_to_log("Chunk {chunk_id} of {service_name} is at second gate".format(
        chunk_id=chunk_id, service_name=service_name))
    with mx:
        if lc[0] is not None:
            sleep_time = (60.0 / MAX_DISPATCHES_PER_MINUTE) - (
                time.time() - lc[0])
            if sleep_time > 0:
                write_to_log(
                    "Sleeping for {:3.1f} seconds to rate-limit wait_for_server_ip.".
                    format(sleep_time))
                time.sleep(sleep_time)
        lc[0] = time.time()
        # if we had to wait here, that counts toward the rate limit delay
        result = wait_for_server_ip_work(service_name, key_path,
                                         remote_username, environment,
                                         max_concurrent, chunk_id)
        return result


# Job functions


def chunk_input(input_files_basenames, chunk_nlines, chunksize):
    """Chunk input files into pieces for performance and parallelism."""
    part_lists = []  # Lists of partial files
    known_nlines = None
    part_suffix = ""

    for input_file in input_files_basenames:
        input_file_full_local_path = result_dir(input_file)

        # Count number of lines in the file
        nlines = int(
            execute_command_with_output(
                "wc -l %s" % input_file_full_local_path).strip().split()[0])
        # Number of lines should be the same in paired files
        if known_nlines is not None:
            msg = "Mismatched line counts in supposedly paired files: {}".format(
                input_files_basenames)
            assert nlines == known_nlines, msg
        known_nlines = nlines

        # Set number of pieces and names
        numparts = (nlines + chunk_nlines - 1) // chunk_nlines
        ndigits = len(str(numparts - 1))
        part_suffix = "-chunksize-%d-numparts-%d-part-" % (chunksize, numparts)
        out_prefix_base = os.path.basename(input_file) + part_suffix
        out_prefix = os.path.join(CHUNKS_RESULT_DIR, out_prefix_base)

        # Split large file into smaller named pieces
        execute_command("split -a %d --numeric-suffixes -l %d %s %s" %
                        (ndigits, chunk_nlines, input_file_full_local_path,
                         out_prefix))
        execute_command(
            "aws s3 cp --quiet %s/ %s/ --recursive --exclude '*' --include '%s*'"
            % (CHUNKS_RESULT_DIR, SAMPLE_S3_OUTPUT_CHUNKS_PATH,
               out_prefix_base))

        # Get the partial file names
        partial_files = []
        paths = execute_command_with_output(
            "ls %s*" % out_prefix).rstrip().split("\n")
        for pf in paths:
            partial_files.append(os.path.basename(pf))

        # Check that the partial files match our expected chunking pattern
        pattern = "{:0%dd}" % ndigits
        expected_partial_files = [(out_prefix_base + pattern.format(i))
                                  for i in range(numparts)]
        msg = "something went wrong with chunking: {} != {}".format(
            partial_files, expected_partial_files)
        assert expected_partial_files == partial_files, msg
        part_lists.append(partial_files)

    # Ex: [["input_R1.fasta-part-1", "input_R2.fasta-part-1"],
    # ["input_R1.fasta-part-2", "input_R2.fasta-part-2"],
    # ["input_R1.fasta-part-3", "input_R2.fasta-part-3"],...]
    input_chunks = [list(part) for part in zip(*part_lists)]
    return part_suffix, input_chunks


def interpret_min_column_number_string(min_column_number_string,
                                       correct_number_of_output_columns,
                                       try_number):
    if min_column_number_string:
        min_column_number = float(min_column_number_string)
        write_to_log(
            "Try no. %d: Smallest number of columns observed in any line was %d"
            % (try_number, min_column_number))
    else:
        write_to_log("Try no. %d: No hits" % try_number)
        min_column_number = correct_number_of_output_columns
    return min_column_number


def add_blank_line(file_path):
    execute_command("echo '' >> %s;" % file_path)


def remove_blank_line(file_path):
    execute_command("sed -i '$ {/^$/d;}' %s" % file_path)


def run_chunk(part_suffix, remote_home_dir, remote_index_dir, remote_work_dir,
              remote_username, input_files, key_path, service, lazy_run):
    """Dispatch a chunk to worker machines for distributed GSNAP or RAPSearch
    group machines and handle their execution.
    """
    assert service in ("gsnap", "rapsearch")

    chunk_id = input_files[0].split(part_suffix)[-1]
    # TODO: Switch to python 3.6 which supports interpolation in string
    # formatting, and we will half the number of lines below.
    srvc = "rapsearch2" if service == "rapsearch" else service
    multihit_basename = "multihit-{service}-out{part_suffix}{chunk_id}.m8".format(
        service=srvc,
        part_suffix=part_suffix,
        chunk_id=chunk_id,
    )
    multihit_local_outfile = os.path.join(CHUNKS_RESULT_DIR, multihit_basename)
    multihit_remote_outfile = os.path.join(remote_work_dir, multihit_basename)
    multihit_s3_outfile = os.path.join(SAMPLE_S3_OUTPUT_CHUNKS_PATH,
                                       multihit_basename)

    base_str = "aws s3 cp --quiet {s3_path}/{input_fa} {remote_work_dir}/{" \
               "input_fa} "
    download_input_from_s3 = " ; ".join(
        base_str.format(
            s3_path=SAMPLE_S3_OUTPUT_CHUNKS_PATH,
            input_fa=input_fa,
            remote_work_dir=remote_work_dir) for input_fa in input_files)

    base_str = "mkdir -p {remote_work_dir} ; {download_input_from_s3} ; "
    if service == "gsnap":
        commands = base_str + "{remote_home_dir}/bin/gsnapl -A m8 --batch=0 --use-shared-memory=0 --gmap-mode=none --npaths=100 --ordered -t 32 --maxsearch=1000 --max-mismatches=40 -D {remote_index_dir} -d nt_k16 {remote_input_files} > {multihit_remote_outfile}"
    else:
        commands = base_str + "/usr/local/bin/rapsearch -d {remote_index_dir}/nr_rapsearch -e -6 -l 10 -a T -b 0 -v 50 -z 24 -q {remote_input_files} -o {multihit_remote_outfile}"

    commands = commands.format(
        remote_work_dir=remote_work_dir,
        download_input_from_s3=download_input_from_s3,
        remote_home_dir=remote_home_dir,
        remote_index_dir=remote_index_dir,
        remote_input_files=" ".join(
            remote_work_dir + "/" + input_fa for input_fa in input_files),
        multihit_remote_outfile=multihit_remote_outfile if service == "gsnap" else multihit_remote_outfile[:-3]
        # Strip the .m8 for RAPSearch as it adds that
    )

    if not lazy_run or not fetch_lazy_result(multihit_s3_outfile,
                                             multihit_local_outfile):
        correct_number_of_output_columns = 12
        min_column_number = 0
        max_tries = 2
        try_number = 1
        instance_ip = ""

        # Check if every row has correct number of columns (12) in the output
        # file on the remote machine
        while min_column_number != correct_number_of_output_columns \
                and try_number <= max_tries:
            write_to_log("waiting for {} server for chunk {}".format(
                service, chunk_id))
            if service == "gsnap":
                max_concurrent = GSNAPL_MAX_CONCURRENT
            else:
                max_concurrent = RAPSEARCH2_MAX_CONCURRENT

            instance_ip = wait_for_server_ip(service, key_path, remote_username,
                                             ENVIRONMENT, max_concurrent, chunk_id)
            write_to_log("starting alignment for chunk %s on %s server %s" %
                         (chunk_id, service, instance_ip))
            execute_command(
                remote_command(commands, key_path, remote_username, instance_ip))

            if service == "gsnap":
                verification_command = "cat %s" % multihit_remote_outfile
            else:
                # For rapsearch, first remove header lines starting with '#'
                verification_command = "grep -v '^#' %s" % multihit_remote_outfile
            verification_command += " | awk '{print NF}' | sort -nu | head -n 1"
            min_column_number_string = execute_command_with_output(
                remote_command(verification_command, key_path, remote_username,
                               instance_ip))
            min_column_number = interpret_min_column_number_string(
                min_column_number_string, correct_number_of_output_columns,
                try_number)
            try_number += 1
            
        # Move output from remote machine to local machine
        msg = "Chunk %s output corrupt; not copying to S3. Re-start pipeline " \
              "to try again." % chunk_id
        assert min_column_number == correct_number_of_output_columns, msg

        with iostream_uploads:  # Limit concurrent uploads so as not to stall the pipeline.
            with iostream:      # Still counts toward the general semaphore.
                execute_command(
                    scp(key_path, remote_username, instance_ip,
                        multihit_remote_outfile, multihit_local_outfile))
                execute_command("aws s3 cp --quiet %s %s/" %
                                (multihit_local_outfile, SAMPLE_S3_OUTPUT_CHUNKS_PATH))
        write_to_log("finished alignment for chunk %s on %s server %s" %
                     (chunk_id, service, instance_ip))
    return multihit_local_outfile


@run_in_subprocess
def call_hits_m8(input_m8, lineage_map_path, accession2taxid_dict_path,
                 output_m8, output_summary):
    """
    Determine the optimal taxon assignment for each read from the alignment
    results. When a read aligns to multiple distinct references, we need to
    assess at which level in the taxonomic hierarchy the multiple alignments
    reach consensus. We refer to this process of controlling for specificity
    as 'hit calling'.

    Input:
    - m8 file of multiple alignments per read

    Outputs:
    - cleaned m8 file with a single, optimal alignment per read
    - file with summary information, including taxonomy level at which
    specificity is reached

    Details:
    - A taxon is a group of any rank (e.g. species, genus, family, etc.).

    - A hit is a match of a read to a known reference labeled with an
    accession ID. We use NCBI's mapping of accession IDs to taxonomy IDs in
    order to retrieve the full taxonomic hierarchy for the accession ID.

    - The full taxonomy hierarchy for a hit is called its "lineage" (species,
    genus, family, etc.). A hit will normally have (positive) NCBI taxon IDs
    at all levels of the hierarchy, but there are some exceptions:

        - We use an artificial negative taxon ID if we have determined that
        the alignment is not specific at the taxonomy level under
        consideration. This happens when a read's multiple reference matches
        do not agree on taxon ID at the given level.

        For example, a read may match 5 references that all belong to
        different species (e.g. Escherichia albertii, Escherichia vulneris,
        Escherichia coli, ...), but to the same genus (Escherichia). In this
        case, we use the taxon ID for the genus (Escherichia) at the
        genus-level, but we populate the species-level with an artificial
        negative ID. The artificial ID is defined based on a negative base (
        INVALID_CALL_BASE_ID), the taxon level (e.g. 2 for genus), and the
        valid parent ID (e.g. genus Escherichia's taxon ID): see helper
        function cleaned_taxid_lineage for the precise formula.

        - Certain entries in NCBI may not have a full lineage classification;
        for example species and family will be defined but genus will be
        undefined. In this case, we populate the undefined taxonomic level
        with an artificial negative ID defined in the same manner as above.
    """
    lineage_map = shelve.open(lineage_map_path)
    accession2taxid_dict = shelve.open(accession2taxid_dict_path)
    # Helper functions
    # TODO: Represent taxids by numbers instead of strings to greatly reduce
    # memory footprint and increase speed.
    lineage_cache = {}

    def get_lineage(accession_id):
        """Find the lineage of the accession ID and utilize a cache for
        performance by reducing random IOPS, ameliorating a key performance
        bottleneck
        """
        if accession_id in lineage_cache:
            return lineage_cache[accession_id]
        accession_taxid = accession2taxid_dict.get(
            accession_id.split(".")[0], "NA")
        result = lineage_map.get(accession_taxid, NULL_LINEAGE)
        lineage_cache[accession_id] = result
        return result

    def accumulate(hits, accession_id):
        """Accumulate hits for summarizing hit information and specificity at
        each taxonomy level
        """
        lineage_taxids = get_lineage(accession_id)
        for level, taxid_at_level in enumerate(lineage_taxids):
            if int(taxid_at_level) < 0:
                # Skip if we have a negative taxid. When an accession doesn't
                # provide species level info, it doesn't contradict any info
                # provided by other accessions. This occurs a lot and
                # handling it in this way seems to work well.
                continue
            hits[level][taxid_at_level] = accession_id

    def call_hit_level(hits):
        for level, hits_at_level in enumerate(hits):
            if len(hits_at_level) == 1:
                taxid, accession_id = hits_at_level.popitem()
                return level + 1, taxid, accession_id
        return -1, "-1", None

    # Read input_m8 and group hits by read id
    m8 = defaultdict(list)
    for read_id, accession_id, _percent_id, _alignment_length, e_value, _line in iterate_m8(
            input_m8, "call_hits_m8_initial_scan"):
        m8[read_id].append((accession_id, e_value))

    # Deduplicate m8 and summarize hits
    summary = {}
    count = 0
    LOG_INCREMENT = 50000
    write_to_log("Starting to summarize hits for {} read ids from {}.".format(
        len(m8), input_m8))
    for read_id, accessions in m8.iteritems():
        # The Expect value (E) is a parameter that describes the number of
        # hits one can 'expect' to see by chance when searching a database of
        # a particular size. It decreases exponentially as the Score (S) of
        # the match increases. Essentially, the E value describes the random
        # background noise. https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web
        # &PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ
        my_best_evalue = min(acc[1] for acc in accessions)
        hits = [{}, {}, {}]
        for accession_id, e_value in accessions:
            if e_value == my_best_evalue:
                accumulate(hits, accession_id)
        summary[read_id] = my_best_evalue, call_hit_level(hits)
        count += 1
        if count % LOG_INCREMENT == 0:
            msg = "Summarized hits for {} read ids from {}, and counting.".format(count, input_m8)
            write_to_log(msg)
    write_to_log("Summarized hits for all {} read ids from {}.".format(
        count, input_m8))

    # Generate output files. outf is the main output_m8 file and outf_sum is
    # the summary level info.
    emitted = set()
    with open(output_m8, "wb") as outf:
        with open(output_summary, "wb") as outf_sum:
            # Iterator over the lines of the m8 file. Emit the hit with the
            # best value that provides the most specific taxonomy
            # information. If there are multiple hits (also called multiple
            # accession IDs) for a given read that all have the same e-value,
            # some may provide species information and some may only provide
            # genus information. We want to emit the one that provides the
            # species information because from that we can infer the rest of
            # the lineage. If we accidentally emitted the one that provided
            # only genus info, downstream steps may have difficulty
            # recovering the species.

            # TODO: Consider all hits within a fixed margin of the best e-value.
            # This change may need to be accompanied by a change to
            # GSNAP/RAPSearch parameters.
            for read_id, accession_id, _percent_id, _alignment_length, e_value, line in iterate_m8(input_m8, "call_hits_m8_emit_deduped_and_summarized_hits"):
                if read_id in emitted:
                    continue

                # Read the fields from the summary level info
                best_e_value, (hit_level, taxid,
                               best_accession_id) = summary[read_id]
                if best_e_value == e_value and best_accession_id in (
                        None, accession_id):
                    # Read out the hit with the best value that provides the
                    # most specific taxonomy information.
                    emitted.add(read_id)
                    outf.write(line)
                    msg = "{read_id}\t{hit_level}\t{taxid}\n".format(
                        read_id=read_id, hit_level=hit_level, taxid=taxid)
                    outf_sum.write(msg)

    # Upload results
    execute_command("aws s3 cp --quiet %s %s/" % (output_m8,
                                                  SAMPLE_S3_OUTPUT_PATH))
    execute_command("aws s3 cp --quiet %s %s/" % (output_summary,
                                                  SAMPLE_S3_OUTPUT_PATH))


def run_remotely(input_files, service, lazy_run):
    """Run a pipeline step / service on remote machines. Set up environment
    for different configurations and chunk files for performance.
    """

    assert service in ("gsnap", "rapsearch")
    key_path = fetch_key(ENVIRONMENT)
    if service == "gsnap":
        remote_username = "ubuntu"
        remote_home_dir = "/home/%s" % remote_username
        remote_work_dir = "%s/batch-pipeline-workdir/%s" % (remote_home_dir,
                                                            SAMPLE_NAME)
        remote_index_dir = "%s/share" % remote_home_dir
        chunk_size = GSNAPL_CHUNK_SIZE
        chunks_in_flight = chunks_in_flight_gsnap
        output_file = MULTIHIT_GSNAPL_OUT
    elif service == "rapsearch":
        remote_username = "ec2-user"
        remote_home_dir = "/home/%s" % remote_username
        remote_work_dir = "%s/data/batch-pipeline-workdir/%s" % (
            remote_home_dir, SAMPLE_NAME)
        remote_index_dir = "%s/references/nr_rapsearch" % remote_home_dir
        chunk_size = RAPSEARCH_CHUNK_SIZE
        chunks_in_flight = chunks_in_flight_rapsearch
        output_file = MULTIHIT_RAPSEARCH_OUT

    # Split files into chunks for performance
    part_suffix, input_chunks = chunk_input(input_files, 2 * chunk_size,
                                            chunk_size)
    # Process chunks
    chunk_output_files = [None] * len(input_chunks)
    chunk_threads = []
    mutex = threading.RLock()
    # Randomize execution order for performance
    randomized = list(enumerate(input_chunks))
    random.shuffle(randomized)

    for n, chunk_input_files in randomized:
        chunks_in_flight.acquire()
        check_for_errors(mutex, chunk_output_files, input_chunks, service)
        t = threading.Thread(
            target=run_chunk_wrapper,
            args=[
                chunks_in_flight, chunk_output_files, n, mutex, run_chunk, [
                    part_suffix, remote_home_dir, remote_index_dir,
                    remote_work_dir, remote_username, chunk_input_files,
                    key_path, service, lazy_run
                ]
            ])
        t.start()
        chunk_threads.append(t)

    # Check chunk completion
    for ct in chunk_threads:
        ct.join()
        check_for_errors(mutex, chunk_output_files, input_chunks, service)
    assert None not in chunk_output_files

    # Concatenate the pieces and upload results
    concatenate_files(chunk_output_files, result_dir(output_file))
    with iostream_uploads:  # Limit concurrent uploads so as not to stall the pipeline.
        with iostream:      # Still counts toward the general semaphore.
            execute_command("aws s3 cp --quiet %s/%s %s/" %
                            (RESULT_DIR, output_file, SAMPLE_S3_OUTPUT_PATH))


def fetch_deuterostome_file(lock=threading.RLock()):  #pylint: disable=dangerous-default-value
    with lock:
        deuterostome_file_basename = os.path.basename(DEUTEROSTOME_TAXIDS)
        deuterostome_file = os.path.join(REF_DIR, deuterostome_file_basename)
        if not os.path.isfile(deuterostome_file):
            execute_command("aws s3 cp --quiet %s %s/" % (DEUTEROSTOME_TAXIDS,
                                                          REF_DIR))
            write_to_log("downloaded deuterostome list")
        return deuterostome_file


def run_chunk_wrapper(chunks_in_flight, chunk_output_files, n, mutex, target,
                      args):
    result = "error"
    try:
        result = target(*args)
    except:
        with mutex:
            traceback.print_exc()
    finally:
        with mutex:
            chunk_output_files[n] = result
        chunks_in_flight.release()


def check_for_errors(mutex, chunk_output_files, input_chunks, task):
    with mutex:
        if "error" in chunk_output_files:
            # We already have per-chunk retries to address transient (non-deterministic) system issues.
            # If a chunk fails on all retries, it must be due to a deterministic & reproducible problem
            # with the chunk data or command options, so we should not even attempt the other chunks.
            err_i = chunk_output_files.index("error")
            raise RuntimeError("All retries failed for {} chunk {}.".format(
                task, input_chunks[err_i]))


def run_generate_unidentified_fasta(input_fa, output_fa):
    #TODO  remove annotated fasta intermediate file and replace > with : below
    subprocess.check_output(
        "grep -A 1 '>NR::NT::' %s | sed '/^--$/d' > %s" % (input_fa,
                                                           output_fa),
        shell=True)
    write_to_log("finished job")
    execute_command("aws s3 cp --quiet %s %s/" % (output_fa,
                                                  SAMPLE_S3_OUTPUT_PATH))


def get_alignment_version_s3_path():
    return os.path.join(SAMPLE_S3_OUTPUT_PATH, VERSION_OUT)


def fetch_input_and_replace_whitespace(input_filename, result):
    s3_input_path = os.path.join(SAMPLE_S3_INPUT_PATH, input_filename)
    s3_output_path = os.path.join(
        SAMPLE_S3_OUTPUT_PATH,
        input_filename)  # usually the same as the input path
    cleaned_input_path = result_dir("nospace.%s" % input_filename)
    try:
        stdouterr = execute_command_with_output(
            "aws s3 cp --quiet {s3_input_path} - | sed 's/[[:blank:]]/{replacement}/g' > {cleaned_input_path}".
            format(
                replacement=";",
                s3_input_path=s3_input_path,
                cleaned_input_path=cleaned_input_path),
            merge_stderr=True)
        # If the file doesn't exist in s3, the above command, sadly, doesn't crash.  But it outputs "download failed" on stderr.
        if "download failed" in stdouterr or os.path.getsize(
                cleaned_input_path) == 0:
            write_to_log(
                "Failed to download {}".format(s3_input_path), warning=True)
            result[0] = None
        else:
            result[0] = cleaned_input_path
    except:
        result[0] = None
        with print_lock:
            traceback.print_exc()
    # This is extremely rare (debug/development only?) and we don't care if it succeeds
    if s3_input_path != s3_output_path:
        threading.Thread(
            target=execute_command,
            args=[
                "aws s3 cp --quiet {s3_input_path} {s3_output_path}".format(
                    s3_input_path=s3_input_path, s3_output_path=s3_output_path)
            ])


def fetch_and_clean_inputs(input_file_list):
    # Fetch inputs and remove tabs in parallel.
    cleaned_inputs = [["ERROR"], ["ERROR"], ["ERROR"]]
    input_fetcher_threads = [
        threading.Thread(
            target=fetch_input_and_replace_whitespace,
            args=[input_file_list[0], cleaned_inputs[0]]),
        threading.Thread(
            target=fetch_input_and_replace_whitespace,
            args=[input_file_list[1], cleaned_inputs[1]]),
        threading.Thread(
            target=fetch_input_and_replace_whitespace,
            args=[input_file_list[2], cleaned_inputs[2]])
    ]
    for ift in input_fetcher_threads:
        ift.start()
    for ift in input_fetcher_threads:
        ift.join()

    # flatten
    cleaned_inputs = [ci[0] for ci in cleaned_inputs]

    for i, ci in enumerate(cleaned_inputs):
        assert ci != "ERROR", "Error fetching input {}".format(i)
        assert ci is None or os.path.isfile(
            ci), "Local file {} not found after download of input {}".format(
                ci, i)

    assert cleaned_inputs[0] is not None, "Input {} not found.  This input is mandatory.".format(
        input_file_list[0])

    assert (cleaned_inputs[1] is None) == (
        cleaned_inputs[2] is
        None), "Input {} is required when {} is given, and vice versa".format(
            input_file_list[1], input_file_list[2])

    cleaned_inputs = filter(None, cleaned_inputs)
    assert len(cleaned_inputs) in (1, 3)

    return cleaned_inputs


def run_stage2(lazy_run=True):
    # Make local directories
    execute_command("mkdir -p %s %s %s %s %s" %
                    (SAMPLE_DIR, FASTQ_DIR, RESULT_DIR, CHUNKS_RESULT_DIR,
                     REF_DIR))

    # Configure logger
    log_file = "%s/%s.%s.txt" % (RESULT_DIR, LOGS_OUT_BASENAME,
                                 AWS_BATCH_JOB_ID)
    configure_logger(log_file)

    write_to_log("Starting stage...")

    # Initiate fetch of references
    lineage_paths = [None]
    accession2taxid_paths = [None]

    def fetch_references():
        lineage_paths[0] = fetch_reference(LINEAGE_SHELF)
        accession2taxid_paths[0] = fetch_reference(ACCESSION2TAXID)

    reference_fetcher_thread = threading.Thread(target=fetch_references)
    reference_fetcher_thread.start()

    # Import existing job stats
    stats = StatsFile(STATS_OUT, RESULT_DIR, SAMPLE_S3_INPUT_PATH,
                      SAMPLE_S3_OUTPUT_PATH)
    stats.load_from_s3()

    if stats.gsnap_ran_in_host_filtering():
        raw_inputs = [
            EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1,
            EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT2,
            EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT3
        ]
    else:
        raw_inputs = [
            EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1,
            EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2,
            EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT3
        ]

    cleaned_inputs = fetch_and_clean_inputs(raw_inputs)

    gsnapl_input_files = [os.path.basename(f) for f in cleaned_inputs[:2]]
    merged_fasta = cleaned_inputs[-1]
    before_file_name_for_log = cleaned_inputs[0]
    before_file_type_for_log = "fasta"  # Unpaired
    if len(gsnapl_input_files) == 2:
        before_file_type_for_log = "fasta_paired"

    # Track if threads succeeded
    thread_success = {}
    thread_success_lock = threading.RLock()
    uploader_threads = {}

    def upload(thread_name, local_path, s3_path):
        execute_command("aws s3 cp --quiet {local_path} {s3_path}".format(
            local_path=local_path, s3_path=s3_path))
        with thread_success_lock:
            thread_success[thread_name] = True

    # Subsample if specified. Use a lock and thread structure to execute the steps with
    # parallelism.
    if SUBSAMPLE:
        target_n_reads = int(SUBSAMPLE)
        if len(gsnapl_input_files) == 2:
            subsampled_gsnapl_input_files, subsampled_merged_fasta = subsample_paired_fastas(
                gsnapl_input_files, os.path.basename(merged_fasta),
                target_n_reads)
        else:
            subsampled_merged_fasta = subsample_single_fasta(
                gsnapl_input_files[0], target_n_reads)
            subsampled_gsnapl_input_files = [subsampled_merged_fasta]

        gsnapl_input_files = subsampled_gsnapl_input_files
        merged_fasta = result_dir(subsampled_merged_fasta)

        for i, f in enumerate(gsnapl_input_files):
            thread_name = "uploader_{}".format(i)
            uploader_threads[thread_name] = threading.Thread(
                target=upload,
                args=[thread_name,
                      result_dir(f), SAMPLE_S3_OUTPUT_PATH + "/"])
            uploader_threads[thread_name].start()

        if len(gsnapl_input_files) == 2:
            thread_name = "uploader_{}".format(len(gsnapl_input_files))
            uploader_threads[thread_name] = threading.Thread(
                target=upload,
                args=[thread_name, merged_fasta, SAMPLE_S3_OUTPUT_PATH + "/"])
            uploader_threads[thread_name].start()

    deuterostome_fetcher = threading.Thread(target=fetch_deuterostome_file)
    deuterostome_fetcher.start()

    def run_gsnap():
        # Run GSNAP remotely
        log_params = return_merged_dict(
            DEFAULT_LOG_PARAMS, {
                "title": "GSNAPL",
                "version_file_s3": GSNAP_VERSION_FILE_S3,
                "output_version_file": result_dir(VERSION_OUT)
            })
        run_and_log_s3(log_params, [result_dir(MULTIHIT_GSNAPL_OUT)], lazy_run,
                       SAMPLE_S3_OUTPUT_PATH, run_remotely, gsnapl_input_files,
                       "gsnap", lazy_run)

        reference_fetcher_thread.join()
        call_hits_m8(
            result_dir(MULTIHIT_GSNAPL_OUT), lineage_paths[0],
            accession2taxid_paths[0], result_dir(DEDUP_MULTIHIT_GSNAPL_OUT),
            result_dir(SUMMARY_MULTIHIT_GSNAPL_OUT))
        stats.count_reads(
            "run_gsnapl_remotely",
            before_filename=before_file_name_for_log,
            before_filetype=before_file_type_for_log,
            after_filename=result_dir(DEDUP_MULTIHIT_GSNAPL_OUT),
            after_filetype="m8")

        deuterostome_path = None
        if not SKIP_DEUTERO_FILTER:
            deuterostome_fetcher.join()
            deuterostome_path = fetch_deuterostome_file()

        generate_taxon_count_json_from_m8(
            result_dir(DEDUP_MULTIHIT_GSNAPL_OUT),
            result_dir(SUMMARY_MULTIHIT_GSNAPL_OUT), 'raw', 'NT',
            lineage_paths[0], deuterostome_path, stats.get_total_reads(),
            stats.get_remaining_reads(), result_dir(MULTIHIT_NT_JSON_OUT))

        with thread_success_lock:
            thread_success["gsnap"] = True

    def run_rapsearch2():
        log_params = return_merged_dict(
            DEFAULT_LOG_PARAMS, {
                "title": "RAPSearch2",
                "version_file_s3": RAPSEARCH_VERSION_FILE_S3,
                "output_version_file": result_dir(VERSION_OUT)
            })
        run_and_log_s3(log_params, [result_dir(MULTIHIT_RAPSEARCH_OUT)],
                       lazy_run, SAMPLE_S3_OUTPUT_PATH, run_remotely,
                       [merged_fasta], "rapsearch", lazy_run)

        reference_fetcher_thread.join()
        call_hits_m8(
            result_dir(MULTIHIT_RAPSEARCH_OUT), lineage_paths[0],
            accession2taxid_paths[0], result_dir(DEDUP_MULTIHIT_RAPSEARCH_OUT),
            result_dir(SUMMARY_MULTIHIT_RAPSEARCH_OUT))
        stats.count_reads(
            "run_rapsearch2_remotely",
            before_filename=result_dir(merged_fasta),
            before_filetype="fasta",
            after_filename=result_dir(DEDUP_MULTIHIT_RAPSEARCH_OUT),
            after_filetype="m8")

        deuterostome_path = None
        if not SKIP_DEUTERO_FILTER:
            deuterostome_fetcher.join()
            deuterostome_path = fetch_deuterostome_file()

        # PRODUCE NEW MULTIHIT NR OUTPUT
        generate_taxon_count_json_from_m8(
            result_dir(DEDUP_MULTIHIT_RAPSEARCH_OUT),
            result_dir(SUMMARY_MULTIHIT_RAPSEARCH_OUT), 'log10', 'NR',
            lineage_paths[0], deuterostome_path, stats.get_total_reads(),
            stats.get_remaining_reads(), result_dir(MULTIHIT_NR_JSON_OUT))

        with thread_success_lock:
            thread_success["rapsearch2"] = True

    def run_additional_steps():

        # COMBINE NEW MULTIHIT NT AND NR OUTPUTS
        combine_pipeline_output_json(
            result_dir(MULTIHIT_NT_JSON_OUT), result_dir(MULTIHIT_NR_JSON_OUT),
            result_dir(MULTIHIT_COMBINED_JSON_OUT), stats)

        with thread_success_lock:
            thread_success["additional_steps"] = True

    def run_annotation_steps():
        # Annotate fasta with accessions
        annotate_fasta_with_accessions(
            merged_fasta, result_dir(DEDUP_MULTIHIT_GSNAPL_OUT),
            result_dir(DEDUP_MULTIHIT_RAPSEARCH_OUT),
            result_dir(ACCESSION_ANNOTATED_FASTA))

        log_params = return_merged_dict(
            DEFAULT_LOG_PARAMS,
            {"title": "generate FASTA of unidentified reads"})
        run_and_log_eager(log_params, [result_dir(UNIDENTIFIED_FASTA_OUT)],
                          run_generate_unidentified_fasta,
                          result_dir(ACCESSION_ANNOTATED_FASTA),
                          result_dir(UNIDENTIFIED_FASTA_OUT))
        stats.count_reads(
            "run_generate_unidentified_fasta",
            before_filename=result_dir(ACCESSION_ANNOTATED_FASTA),
            before_filetype="fasta",
            after_filename=result_dir(UNIDENTIFIED_FASTA_OUT),
            after_filetype="fasta")

        with thread_success_lock:
            thread_success["annotation"] = True

    t_gsnap = threading.Thread(target=run_gsnap)
    t_rapsearch2 = threading.Thread(target=run_rapsearch2)
    t_additional_steps = threading.Thread(target=run_additional_steps)
    t_annotation = threading.Thread(target=run_annotation_steps)

    def thread_succeeded(name):
        with thread_success_lock:
            return thread_success.get(name)

    # Main steps
    t_gsnap.start()
    t_rapsearch2.start()
    t_gsnap.join()
    assert thread_succeeded("gsnap")
    t_rapsearch2.join()
    assert thread_succeeded("rapsearch2")

    # Additional steps
    t_additional_steps.start()
    t_annotation.start()
    t_additional_steps.join()
    assert thread_succeeded("additional_steps")
    t_annotation.join()
    assert thread_succeeded("annotation")

    for thread_name, t in uploader_threads.items():
        t.join()
        assert thread_succeeded(thread_name), "thread {} failed".format(
            thread_name)

    # copy log file -- after work is done
    stats.save_to_s3()
    upload_log_file(SAMPLE_S3_OUTPUT_PATH)
