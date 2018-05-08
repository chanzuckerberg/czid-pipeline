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
from .common import * #pylint: disable=wildcard-import


# that's per job;  by default in early 2018 jobs are subsampled to <= 100 chunks
MAX_CHUNKS_IN_FLIGHT_GSNAP = 32
MAX_CHUNKS_IN_FLIGHT_RAPSEARCH = 32
chunks_in_flight_gsnap = threading.Semaphore(MAX_CHUNKS_IN_FLIGHT_GSNAP)
chunks_in_flight_rapsearch = threading.Semaphore(MAX_CHUNKS_IN_FLIGHT_RAPSEARCH)

# Dispatch at most this many chunks per minute, to ensure fairness
# amongst jobs regardless of job size (not the best way to do it,
# just for now).
MAX_DISPATCHES_PER_MINUTE = 10

# poll this many random servers per wait_for_instance_ip
MAX_INSTANCES_TO_POLL = 8

MAX_POLLING_LATENCY = 10 # seconds

# If no instance is available, we should refresh our list of instances, to pick up new instances
# added by autoscaling.  Wait at least this long between refreshes to stay within AWS account limits.
MIN_INTERVAL_BETWEEN_DESCRIBE_INSTANCES = 180

# Refresh at least every 30 minutes
MAX_INTERVAL_BETWEEN_DESCRIBE_INSTANCES = 900

# data directories
ROOT_DIR = '/mnt'
DEST_DIR = ROOT_DIR + '/idseq/data' # generated data go here
REF_DIR = ROOT_DIR + '/idseq/ref' # referene genome / ref databases go here
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
SUBSAMPLE = os.environ.get('SUBSAMPLE') # number of read pairs to subsample to, before gsnap/rapsearch
FASTQ_BUCKET = os.environ.get('FASTQ_BUCKET')
INPUT_BUCKET = os.environ.get('INPUT_BUCKET')
FILE_TYPE = os.environ.get('FILE_TYPE', 'fastq.gz')
OUTPUT_BUCKET = os.environ.get('OUTPUT_BUCKET')
AWS_BATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID', 'local')
ENVIRONMENT = os.environ.get('ENVIRONMENT', 'production')
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
DEFAULT_LOGPARAMS = {"sample_s3_output_path": SAMPLE_S3_OUTPUT_PATH}

# For reproducibility of random operations
random.seed(hash(SAMPLE_NAME))

# versioning
## For now, index updates are infrequent and we can get their versions from S3.
## If updates ever become frequent, we may want to check instead which version is actually on the
## machine taking the job, possibly out of sync with the newest version in S3.

base_s3 = 's3://idseq-database/alignment_indexes'
base_dt = '2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime'
GSNAP_VERSION_FILE_S3 = ("%s/%s/nt_k16.version.txt" % (base_s3, base_dt))
RAPSEARCH_VERSION_FILE_S3 = ("%s/%s/nr_rapsearch.version.txt" % (base_s3, base_dt))

# compute capacity
GSNAPL_MAX_CONCURRENT = 3 # number of gsnapl jobs allowed to run concurrently on 1 machine
RAPSEARCH2_MAX_CONCURRENT = 6
GSNAPL_CHUNK_SIZE = 15000 # number of fasta records in a chunk so it runs in ~10 minutes on i3.16xlarge
RAPSEARCH_CHUNK_SIZE = 10000

# references
# from common import ACCESSION2TAXID
base_s3 = 's3://idseq-database/taxonomy'
base_dt = '2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime'
DEUTEROSTOME_TAXIDS = ("%s/%s/deuterostome_taxids.txt" % (base_s3, base_dt))
# from common import LINEAGE_SHELF

# definitions for integration with web app
TAX_LEVEL_SPECIES = 1
TAX_LEVEL_GENUS = 2
TAX_LEVEL_FAMILY = 3
MISSING_GENUS_ID = -200
MISSING_FAMILY_ID = -300


# convenience functions
def count_lines(input_file):
    return int(execute_command_with_output("wc -l %s" % input_file).strip().split()[0])


def count_lines_in_paired_files(input_files):
    distinct_counts = list(set(map(count_lines, input_files)))
    assert len(distinct_counts) == 1, "Mismatched line counts in supposedly paired files: {}".format(input_files)
    return distinct_counts[0]


def subsample_helper(input_file, records_to_keep, type_, output_file):
    record_number = 0
    assert isinstance(records_to_keep, set), "Essential for this to complete in our lifetimes."
    kept_read_ids = set()
    kept_reads_count = 0
    # write_to_log(sorted(records_to_keep))
    with open(input_file, 'rb') as input_f:
        with open(output_file, 'wb') as output_f:
            sequence_name = input_f.readline()
            sequence_data = input_f.readline()
            while len(sequence_name) > 0 and len(sequence_data) > 0:
                if type_ == "read_ids":
                    sequence_basename = sequence_name.rstrip().rsplit('/', 1)[0]
                    condition = sequence_basename in records_to_keep
                elif type_ == "record_indices":
                    condition = record_number in records_to_keep
                    sequence_basename = sequence_name.rstrip()
                else:
                    condition = True
                if condition:
                    output_f.write(sequence_name)
                    output_f.write(sequence_data)
                    kept_reads_count += 1
                    kept_read_ids.add(sequence_basename)
                sequence_name = input_f.readline()
                sequence_data = input_f.readline()
                record_number += 1
    if type_ == "record_indices":
        if len(kept_read_ids) != len(records_to_keep):
            write_to_log("WARNING:  Repeated read IDs in {input_file}:  Found {kri} read IDs in {rtk} sampled rows.".format(input_file=input_file, kri=len(kept_read_ids), rtk=len(records_to_keep)))
    if type_ == "read_ids":
        if kept_read_ids != records_to_keep:
            missing = records_to_keep - kept_read_ids
            examples = sorted(random.sample(missing, min(10, len(missing))))
            assert kept_read_ids == records_to_keep, "Not all desired read IDs were found in the file: {}\nMissing: {}".format(input_file, examples)
    return kept_read_ids, kept_reads_count


# Note dedicated random stream for just this function, so that arbitrary other use of randomness (e.g. for I/O retry delays) will not perturb the subsampling stream
def subsample_single_fasta(input_files_basename, target_n_reads, randgen=random.Random(x=hash(SAMPLE_NAME))):
    input_file = result_dir(input_files_basename)
    total_records = count_lines(input_file) // 2 # each fasta record spans 2 lines
    write_to_log("total unpaired reads: %d" % total_records)
    write_to_log("target reads: %d" % target_n_reads)
    if total_records <= target_n_reads:
        return input_file
    records_to_keep = set(randgen.sample(xrange(total_records), target_n_reads))
    subsample_prefix = "subsample_%d" % target_n_reads
    input_dir, input_basename = os.path.split(input_file)
    output_basename = "%s.%s" % (subsample_prefix, input_basename)
    output_file = os.path.join(input_dir, output_basename)
    subsample_helper(input_file, records_to_keep, "record_indices", output_file)
    return output_file


def subsample_paired_fastas(input_files_basenames, merged_file_basename, target_n_reads, randgen=random.Random(x=hash(SAMPLE_NAME))):
    input_files = [result_dir(f) for f in input_files_basenames]
    merged_file = result_dir(merged_file_basename)
    total_records = count_lines_in_paired_files(input_files) // 2 # each fasta record spans 2 lines
    write_to_log("total read pairs: %d" % total_records)
    write_to_log("target read pairs: %d" % target_n_reads)
    # note: target_n_reads and total_records really refer to numbers of read PAIRS
    if total_records <= target_n_reads:
        return input_files_basenames, merged_file_basename
    subsample_prefix = "subsample_%d" % target_n_reads
    records_to_keep = set(randgen.sample(xrange(total_records), target_n_reads))
    subsampled_files = []
    known_kept_read_ids = set()
    kept_reads_count = 0
    # subsample the paired files and record read IDs kept
    for input_file in input_files:
        input_dir = os.path.split(input_file)[0]
        input_basename = os.path.split(input_file)[1]
        output_basename = "%s.%s" % (subsample_prefix, input_basename)
        output_file = os.path.join(input_dir, output_basename)
        kept_read_ids, kept_reads_count = subsample_helper(input_file, records_to_keep, "record_indices", output_file)
        assert kept_reads_count == target_n_reads
        subsampled_files.append(output_basename)
        known_kept_read_ids |= kept_read_ids
    # subsample the merged file to the same read IDs
    input_dir = os.path.split(merged_file)[0]
    input_basename = os.path.split(merged_file)[1]
    subsampled_merged_file = "%s.%s" % (subsample_prefix, input_basename)
    output_file = os.path.join(input_dir, subsampled_merged_file)
    _, kept_reads_count = subsample_helper(merged_file, known_kept_read_ids, "read_ids", output_file)
    num_desired = target_n_reads * len(input_files_basenames)
    if kept_reads_count != num_desired:
        # TODO  Surface these warnings in the webapp.
        write_to_log("ERROR:  Improperly paired reads.  Total reads in sampled merged fasta: {krc}".format(krc=kept_reads_count))
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
    invalid_hits = 0
    last_invalid_line = None
    with open(m8_file, 'rb') as m8f:
        line_count = 0
        for line in m8f:
            line_count += 1
            if line and line[0] == '#':
                # skip comments
                # write_to_log("{}:{}: {}".format(m8_file, line_count, line.rstrip()), flush=False)
                continue
            parts = line.split("\t")
            assert len(parts) >= 12, log_corrupt(m8_file, line)
            read_id = parts[0]
            accession_id = parts[1]
            percent_id = float(parts[2])
            alignment_length = int(parts[3])
            e_value = float(parts[10])
            # GSNAP outputs these sometimes, and usually they are not the only assingment, so rather
            # than killing the job, we just skip them.  If we don't filter these out here, they will
            # override the good data when computing min(evalue), pollute averages computed in the json,
            # and cause the webapp loader to crash as the RAILS JSON parser cannot handle NaNs.
            if alignment_length <= 0 or not -0.25 < percent_id < 100.25 or e_value != e_value:
                invalid_hits += 1
                last_invalid_line = line
                continue
            if debug_caller and line_count % logging_interval == 0:
                write_to_log("Scanned {} m8 lines from {} for {}, and going.".format(line_count, m8_file, debug_caller))
            yield (read_id, accession_id, percent_id, alignment_length, e_value, line)
    if invalid_hits:
        write_to_log("Found {} invalid hits in {};  last invalid hit line: {}".format(invalid_hits, m8_file, last_invalid_line), warning=True)
    if debug_caller:
        write_to_log("Scanned all {} m8 lines from {} for {}.".format(line_count, m8_file, debug_caller))


def annotate_fasta_with_accessions(merged_input_fasta, nt_m8, nr_m8, output_fasta):

    def get_map(m8_file):
        return dict((read_id, accession_id) for read_id, accession_id, _percent_id, _alignment_length, _e_value, _line in iterate_m8(m8_file, "annotate_fasta_with_accessions"))

    nt_map = get_map(nt_m8)
    nr_map = get_map(nr_m8)

    with open(merged_input_fasta, 'rb') as input_fasta_f:
        with open(output_fasta, 'wb') as output_fasta_f:
            sequence_name = input_fasta_f.readline()
            sequence_data = input_fasta_f.readline()
            while sequence_name and sequence_data:
                read_id = sequence_name.rstrip().lstrip('>')
                # need to annotate NR then NT in this order for alignment viz
                new_read_name = "NR:{nr_accession}:NT:{nt_accession}:{read_id}".format(
                    nr_accession=nr_map.get(read_id, ''),
                    nt_accession=nt_map.get(read_id, ''),
                    read_id=read_id
                )
                output_fasta_f.write(">%s\n" % new_read_name)
                output_fasta_f.write(sequence_data)
                sequence_name = input_fasta_f.readline()
                sequence_data = input_fasta_f.readline()
    async_handler.launch_aws_upload(output_fasta, SAMPLE_S3_OUTPUT_PATH + "/")

@run_in_subprocess
def generate_taxon_count_json_from_m8(m8_file, hit_level_file, e_value_type, count_type, lineage_map_path, deuterostome_path, total_reads, remaining_reads, output_file):
    if SKIP_DEUTERO_FILTER:
        def any_hits_to_remove(_hits):
            return False
    else:
        taxids_to_remove = read_file_into_set(deuterostome_path)
        def any_hits_to_remove(hits):
            for taxid in hits:
                if int(taxid) >= 0 and taxid in taxids_to_remove:
                    return True
            return False
    aggregation = {}
    hit_f = open(hit_level_file, 'rb')
    m8_f = open(m8_file, 'rb')
    # lines in m8_file and hit_level_file correspond (same read_id)
    hit_line = hit_f.readline()
    m8_line = m8_f.readline()
    lineage_map = shelve.open(lineage_map_path)
    LINEAGE_WILDCARDS = ("-100", "-200", "-300")
    NUM_RANKS = len(LINEAGE_WILDCARDS)
    # See https://en.wikipedia.org/wiki/Double-precision_floating-point_format
    MIN_NORMAL_POSITIVE_DOUBLE = 2.0 ** -1022
    while hit_line and m8_line:
        # Retrieve data values from files:
        hit_line_columns = hit_line.rstrip("\n").split("\t")
        _read_id = hit_line_columns[0]
        hit_level = hit_line_columns[1]
        hit_taxid = hit_line_columns[2]
        if int(hit_level) < 0:
            hit_line = hit_f.readline()
            m8_line = m8_f.readline()
            continue
        m8_line_columns = m8_line.split("\t")
        assert m8_line_columns[0] == hit_line_columns[0], "read_ids in %s and %s do not match: %s vs. %s" % (os.path.basename(m8_file), os.path.basename(hit_level_file), m8_line_columns[0], hit_line_columns[0])
        percent_identity = float(m8_line_columns[2])
        alignment_length = float(m8_line_columns[3])
        e_value = float(m8_line_columns[10])
        # These have been filtered out before the creation of m8_f and hit_f
        assert alignment_length > 0
        assert -0.25 < percent_identity < 100.25
        assert e_value == e_value
        if e_value_type != 'log10':
            assert MIN_NORMAL_POSITIVE_DOUBLE <= e_value  # Subtle, but this excludes positive subnorms.
            e_value = math.log10(e_value)

        # Retrieve the taxon lineage and mark meaningless calls with fake taxids.
        hit_taxids_all_levels = lineage_map.get(hit_taxid, LINEAGE_WILDCARDS)
        cleaned_hit_taxids_all_levels = validate_taxid_lineage(hit_taxids_all_levels, hit_taxid, hit_level)
        assert NUM_RANKS == len(cleaned_hit_taxids_all_levels)

        if not any_hits_to_remove(cleaned_hit_taxids_all_levels):
            # Aggregate each level
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
        taxon_counts_attributes.append({# TODO: Extend taxonomic ranks as indicated on the commented out lines.
            "tax_id": agg_key[0],
            "tax_level": tax_level,
            # 'species_taxid' : agg_key[tax_level - 1] if tax_level == 1 else "-100",
            'genus_taxid': agg_key[2 - tax_level] if tax_level <= 2 else "-200",
            'family_taxid': agg_key[3 - tax_level] if tax_level <= 3 else "-300",
            # 'order_taxid' : agg_key[4 - tax_level] if tax_level <= 4 else "-400",
            # 'class_taxid' : agg_key[5 - tax_level] if tax_level <= 5 else "-500",
            # 'phyllum_taxid' : agg_key[6 - tax_level] if tax_level <= 6 else "-600",
            # 'kingdom_taxid' : agg_key[7 - tax_level] if tax_level <= 7 else "-700",
            # 'domain_taxid' : agg_key[8 - tax_level] if tax_level <= 8 else "-800",
            "count": count,
            "percent_identity": agg_bucket['sum_percent_identity'] / count,
            "alignment_length": agg_bucket['sum_alignment_length'] / count,
            "e_value": agg_bucket['sum_e_value'] / count,
            "count_type": count_type})
    output_dict = {
        "pipeline_output": {
            "total_reads": total_reads,
            "remaining_reads": remaining_reads,
            "taxon_counts_attributes": taxon_counts_attributes
        }
    }
    with open(output_file, 'wb') as outf:
        json.dump(output_dict, outf)
    async_handler.launch_aws_upload(output_file, SAMPLE_S3_OUTPUT_PATH + "/")


def combine_pipeline_output_json(inputPath1, inputPath2, outputPath, stats):
    total_reads = stats.get_total_reads()
    remaining_reads = stats.get_remaining_reads()
    with open(inputPath1) as inf1:
        input1 = json.load(inf1).get("pipeline_output")
    with open(inputPath2) as inf2:
        input2 = json.load(inf2).get("pipeline_output")
    taxon_counts_attributes = (input1.get("taxon_counts_attributes")
                               + input2.get("taxon_counts_attributes"))
    pipeline_output_dict = {
        "total_reads": total_reads,
        "remaining_reads": remaining_reads,
        "taxon_counts_attributes": taxon_counts_attributes
    }
    output_dict = {
        "pipeline_output": pipeline_output_dict
    }
    with open(outputPath, 'wb') as outf:
        json.dump(output_dict, outf)
    async_handler.launch_aws_upload(outputPath, SAMPLE_S3_OUTPUT_PATH + "/")


def read_file_into_set(file_name):
    with open(file_name) as f:
        S = set(x.rstrip() for x in f)
    S.discard('')
    return S


def environment_for_aligners(_environment):
    return "production" ## temporary fix since we only have "production" gsnap/rapsearch machines right now


def fetch_key(environment, mutex=threading.RLock()):
    with mutex:
        key_s3_path = "s3://idseq-secrets/idseq-%s.pem" % environment_for_aligners(environment)
        key_name = os.path.basename(key_s3_path)
        key_path = REF_DIR +'/' + key_name
        if not os.path.exists(key_path):
            execute_command("aws s3 cp --quiet %s %s/" % (key_s3_path, REF_DIR))
            execute_command("chmod 400 %s" % key_path)
        return key_path


@retry
def get_server_ips_work(service_name, environment):
    tag = "service"
    value = "%s-%s" % (service_name, environment_for_aligners(environment))
    describe_json = json.loads(execute_command_with_output("aws ec2 describe-instances --filters 'Name=tag:%s,Values=%s' 'Name=instance-state-name,Values=running'" % (tag, value)))
    server_ips = []
    for reservation in describe_json["Reservations"]:
        for instance in reservation["Instances"]:
            server_ips += [instance["NetworkInterfaces"][0]["Association"]["PublicIp"]]
    return server_ips


def get_server_ips(service_name, environment, aggressive=False, cache={}, mutex=threading.RLock()): #pylint: disable=dangerous-default-value
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
                cache[cache_key] = (now, get_server_ips_work(service_name, environment))
            return cache[cache_key][1]
    except:
        # return [] causes a sleep of wait_seconds before retrying (see below)
        traceback.print_exc()
        return []


def wait_for_server_ip_work(service_name, key_path, remote_username, environment, max_concurrent, chunk_id, had_to_wait=[False]): #pylint: disable=dangerous-default-value
    while True:
        write_to_log("Chunk {chunk_id} of {service_name} is at third gate".format(chunk_id=chunk_id, service_name=service_name))
        instance_ips = get_server_ips(service_name, environment, aggressive=had_to_wait[0])
        instance_ips = random.sample(instance_ips, min(MAX_INSTANCES_TO_POLL, len(instance_ips)))
        ip_nproc_dict = {}
        dict_mutex = threading.RLock()
        dict_writable = True
        def poll_server(ip):
            command = 'ssh -o "StrictHostKeyChecking no" -o "ConnectTimeout 5" -i %s %s@%s "ps aux | grep %s | grep -v bash" || echo "error"' % (key_path, remote_username, ip, service_name)
            output = execute_command_with_output(command, timeout=MAX_POLLING_LATENCY).rstrip().split("\n")
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
            write_to_log("Only {} out of {} instances responded to polling;  {} threads are timing out.".format(len(ip_nproc_dict), len(poller_threads), len([pt for pt in poller_threads if pt.isAlive()])))
        write_to_log("Chunk {chunk_id} of {service_name} is at fourth gate".format(chunk_id=chunk_id, service_name=service_name))
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
            write_to_log("%s server %s has capacity %d. Kicking off " % (service_name, min_nproc_ip, free_slots))
            return min_nproc_ip
        else:
            had_to_wait[0] = True
            wait_seconds = random.randint(max(20, MAX_POLLING_LATENCY), max(60, MAX_POLLING_LATENCY))
            write_to_log("%s servers busy. Wait for %d seconds" % (service_name, wait_seconds))
            time.sleep(wait_seconds)


def wait_for_server_ip(service_name, key_path, remote_username, environment, max_concurrent, chunk_id, mutex=threading.RLock(), mutexes={}, last_checks={}): #pylint: disable=dangerous-default-value
    # We rate limit these to ensure fairness across jobs regardless of job size
    with mutex:
        if service_name not in mutexes:
            mutexes[service_name] = threading.RLock()
            last_checks[service_name] = [None]
        lc = last_checks[service_name]
        mx = mutexes[service_name]
    write_to_log("Chunk {chunk_id} of {service_name} is at second gate".format(chunk_id=chunk_id, service_name=service_name))
    with mx:
        if lc[0] != None:
            sleep_time = (60.0 / MAX_DISPATCHES_PER_MINUTE) - (time.time() - lc[0])
            if sleep_time > 0:
                write_to_log("Sleeping for {:3.1f} seconds to rate-limit wait_for_server_ip.".format(sleep_time))
                time.sleep(sleep_time)
        lc[0] = time.time()
        # if we had to wait here, that counts toward the rate limit delay
        result = wait_for_server_ip_work(service_name, key_path, remote_username, environment, max_concurrent, chunk_id)
        return result


# job functions
def chunk_input(input_files_basenames, chunk_nlines, chunksize):
    part_lists = []
    known_nlines = None
    for input_file in input_files_basenames:
        input_file_full_local_path = result_dir(input_file)
        nlines = int(execute_command_with_output("wc -l %s" % input_file_full_local_path).strip().split()[0])
        if known_nlines != None:
            assert nlines == known_nlines, "Mismatched line counts in supposedly paired files: {}".format(input_files_basenames)
        known_nlines = nlines
        numparts = (nlines + chunk_nlines - 1) // chunk_nlines
        ndigits = len(str(numparts - 1))
        part_suffix = "-chunksize-%d-numparts-%d-part-" % (chunksize, numparts)
        out_prefix_base = os.path.basename(input_file) + part_suffix
        out_prefix = os.path.join(CHUNKS_RESULT_DIR, out_prefix_base)
        execute_command("split -a %d --numeric-suffixes -l %d %s %s" % (ndigits, chunk_nlines, input_file_full_local_path, out_prefix))
        execute_command("aws s3 cp --quiet %s/ %s/ --recursive --exclude '*' --include '%s*'" % (CHUNKS_RESULT_DIR, SAMPLE_S3_OUTPUT_CHUNKS_PATH, out_prefix_base))
        # note: no sleep(10) after this upload to s3... not sure if that was ever needed
        partial_files = [os.path.basename(partial_file) for partial_file in execute_command_with_output("ls %s*" % out_prefix).rstrip().split("\n")]
        pattern = "{:0%dd}" % ndigits
        expected_partial_files = [(out_prefix_base + pattern.format(i)) for i in range(numparts)]
        assert expected_partial_files == partial_files, "something went wrong with chunking: {} != {}".format(partial_files, expected_partial_files)
        part_lists.append(partial_files)
    input_chunks = [list(part) for part in zip(*part_lists)]
    # e.g. [["input_R1.fasta-part-1", "input_R2.fasta-part-1"],["input_R1.fasta-part-2", "input_R2.fasta-part-2"],["input_R1.fasta-part-3", "input_R2.fasta-part-3"],...]
    return part_suffix, input_chunks


def interpret_min_column_number_string(min_column_number_string, correct_number_of_output_columns, try_number):
    if min_column_number_string:
        min_column_number = float(min_column_number_string)
        write_to_log("Try no. %d: Smallest number of columns observed in any line was %d" % (try_number, min_column_number))
    else:
        write_to_log("Try no. %d: No hits" % try_number)
        min_column_number = correct_number_of_output_columns
    return min_column_number


def add_blank_line(file_path):
    execute_command("echo '' >> %s;" % file_path)


def remove_blank_line(file_path):
    execute_command("sed -i '$ {/^$/d;}' %s" % file_path)


def run_chunk(part_suffix, remote_home_dir, remote_index_dir, remote_work_dir, remote_username, input_files, key_path, service, lazy_run):
    assert service in ("gsnap", "rapsearch")
    chunk_id = input_files[0].split(part_suffix)[-1]
    # TODO: Switch to python 3.6 which supports interpolation in string formatting, and we will half the number of lines below.
    multihit_basename = "multihit-{service}-out{part_suffix}{chunk_id}.m8".format(
        service="rapsearch2" if service == "rapsearch" else service,
        part_suffix=part_suffix,
        chunk_id=chunk_id,
    )
    multihit_local_outfile = os.path.join(CHUNKS_RESULT_DIR, multihit_basename)
    multihit_remote_outfile = os.path.join(remote_work_dir, multihit_basename)
    multihit_s3_outfile = os.path.join(SAMPLE_S3_OUTPUT_CHUNKS_PATH, multihit_basename)
    download_input_from_s3 = " ; ".join(
        "aws s3 cp --quiet {s3_path}/{input_fa} {remote_work_dir}/{input_fa}".format(
            s3_path=SAMPLE_S3_OUTPUT_CHUNKS_PATH,
            input_fa=input_fa,
            remote_work_dir=remote_work_dir)
        for input_fa in input_files
    )
    if service == "gsnap":
        commands = "mkdir -p {remote_work_dir} ; \
            {download_input_from_s3} ; \
            {remote_home_dir}/bin/gsnapl -A m8 --batch=0 --use-shared-memory=0 --gmap-mode=none --npaths=100 --ordered -t 32 --maxsearch=1000 --max-mismatches=40 -D {remote_index_dir} -d nt_k16 {remote_input_files} > {multihit_remote_outfile}"
    else:
        commands = "mkdir -p {remote_work_dir} ; \
            {download_input_from_s3} ; \
            /usr/local/bin/rapsearch -d {remote_index_dir}/nr_rapsearch -e -6 -l 10 -a T -b 0 -v 50 -z 24 -q {remote_input_files} -o {multihit_remote_outfile}"
    commands = commands.format(
        remote_work_dir=remote_work_dir,
        download_input_from_s3=download_input_from_s3,
        remote_home_dir=remote_home_dir,
        remote_index_dir=remote_index_dir,
        remote_input_files=" ".join(remote_work_dir + "/" + input_fa for input_fa in input_files),
        multihit_remote_outfile=multihit_remote_outfile if service == "gsnap" else multihit_remote_outfile[:-3]  # strip the .m8 for rapsearch as it adds that
    )
    if not lazy_run or not fetch_lazy_result(multihit_s3_outfile, multihit_local_outfile):
        correct_number_of_output_columns = 12
        min_column_number = 0
        max_tries = 2
        try_number = 1
        while min_column_number != correct_number_of_output_columns and try_number <= max_tries:
            write_to_log("waiting for {} server for chunk {}".format(service, chunk_id))
            max_concurrent = GSNAPL_MAX_CONCURRENT if service == "gsnap" else RAPSEARCH2_MAX_CONCURRENT
            instance_ip = wait_for_server_ip(service, key_path, remote_username, ENVIRONMENT, max_concurrent, chunk_id)
            write_to_log("starting alignment for chunk %s on %s server %s" % (chunk_id, service, instance_ip))
            execute_command(remote_command(commands, key_path, remote_username, instance_ip))
            # check if every row has correct number of columns (12) in the output file on the remote machine
            if service == "gsnap":
                verification_command = "cat %s" % multihit_remote_outfile
            else:
                # for rapsearch, first, remove header lines starting with '#'
                verification_command = "grep -v '^#' %s" % multihit_remote_outfile
            verification_command += " | awk '{print NF}' | sort -nu | head -n 1"
            min_column_number_string = execute_command_with_output(remote_command(verification_command, key_path, remote_username, instance_ip))
            min_column_number = interpret_min_column_number_string(min_column_number_string, correct_number_of_output_columns, try_number)
            try_number += 1
        # move output from remote machine to local
        assert min_column_number == correct_number_of_output_columns, "Chunk %s output corrupt; not copying to S3. Re-start pipeline to try again." % chunk_id
        execute_command(scp(key_path, remote_username, instance_ip, multihit_remote_outfile, multihit_local_outfile))
        async_handler.launch_aws_upload(multihit_local_outfile, SAMPLE_S3_OUTPUT_CHUNKS_PATH + "/")
        write_to_log("finished alignment for chunk %s on %s server %s" % (chunk_id, service, instance_ip))
    return multihit_local_outfile


@run_in_subprocess
def call_hits_m8(input_m8, lineage_map_path, accession2taxid_dict_path, output_m8, output_summary):
    lineage_map = shelve.open(lineage_map_path)
    accession2taxid_dict = shelve.open(accession2taxid_dict_path)
    # Helper functions
    # TODO: Represent taxids by numbers instead of strings to greatly reduce memory footprint
    # and increase speed.
    NULL_TAXIDS = ("-100", "-200", "-300")
    lineage_cache = {}
    def get_lineage(accession_id):
        if accession_id in lineage_cache:
            result = lineage_cache[accession_id]
        else:
            accession_taxid = accession2taxid_dict.get(accession_id.split(".")[0], "NA")
            result = lineage_map.get(accession_taxid, NULL_TAXIDS)
            lineage_cache[accession_id] = result
        return result
    def accumulate(hits, accession_id):
        lineage_taxids = get_lineage(accession_id)
        for level, taxid_at_level in enumerate(lineage_taxids):
            if int(taxid_at_level) < 0:
                # When an accession doesn't provide species level info, it doesn't contradict
                # any info provided by other accessions.  This occurs a lot and handling it
                # in this way seems to work well.
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
    for read_id, accession_id, _percent_id, _alignment_length, e_value, _line in iterate_m8(input_m8, "call_hits_m8_initial_scan"):
        m8[read_id].append((accession_id, e_value))
    # Deduplicate m8 and summarize hits
    summary = {}
    count = 0
    LOG_INCREMENT = 50000
    write_to_log("Starting to summarize hits for {} read ids from {}.".format(len(m8), input_m8))
    for read_id, accessions in m8.iteritems():
        my_best_evalue = min(acc[1] for acc in accessions)
        hits = [{}, {}, {}]
        for accession_id, e_value in accessions:
            if e_value == my_best_evalue:
                accumulate(hits, accession_id)
        summary[read_id] = my_best_evalue, call_hit_level(hits)
        count += 1
        if count % LOG_INCREMENT == 0:
            write_to_log("Summarized hits for {} read ids from {}, and counting.".format(count, input_m8))
    write_to_log("Summarized hits for all {} read ids from {}.".format(count, input_m8))
    emitted = set()
    with open(output_m8, "wb") as outf:
        with open(output_summary, "wb") as outf_sum:
            for read_id, accession_id, _percent_id, _alignment_length, e_value, line in iterate_m8(input_m8, "call_hits_m8_emit_deduped_and_summarized_hits"):
                if read_id not in emitted:
                    best_e_value, (hit_level, taxid, best_accession_id) = summary[read_id]
                    if best_e_value == e_value and best_accession_id in (None, accession_id):
                        emitted.add(read_id)
                        outf.write(line)
                        outf_sum.write("{read_id}\t{hit_level}\t{taxid}\n".format(read_id=read_id, hit_level=hit_level, taxid=taxid))
    async_handler.launch_aws_upload(output_m8, SAMPLE_S3_OUTPUT_PATH + "/")
    async_handler.launch_aws_upload(output_summary, SAMPLE_S3_OUTPUT_PATH + "/")


def run_remotely(input_files, service, lazy_run):
    assert service in ("gsnap", "rapsearch")
    key_path = fetch_key(ENVIRONMENT)
    if service == "gsnap":
        remote_username = "ubuntu"
        remote_home_dir = "/home/%s" % remote_username
        remote_work_dir = "%s/batch-pipeline-workdir/%s" % (remote_home_dir, SAMPLE_NAME)
        remote_index_dir = "%s/share" % remote_home_dir
        chunk_size = GSNAPL_CHUNK_SIZE
        chunks_in_flight = chunks_in_flight_gsnap
        output_file = MULTIHIT_GSNAPL_OUT
    elif service == "rapsearch":
        remote_username = "ec2-user"
        remote_home_dir = "/home/%s" % remote_username
        remote_work_dir = "%s/data/batch-pipeline-workdir/%s" % (remote_home_dir, SAMPLE_NAME)
        remote_index_dir = "%s/references/nr_rapsearch" % remote_home_dir
        chunk_size = RAPSEARCH_CHUNK_SIZE
        chunks_in_flight = chunks_in_flight_rapsearch
        output_file = MULTIHIT_RAPSEARCH_OUT
    # split file:
    part_suffix, input_chunks = chunk_input(input_files, 2 * chunk_size, chunk_size)
    # process chunks:
    chunk_output_files = [None] * len(input_chunks)
    chunk_threads = []
    mutex = threading.RLock()
    iii = list(enumerate(input_chunks))
    random.shuffle(iii)
    for n, chunk_input_files in iii:
        chunks_in_flight.acquire()
        check_for_errors(mutex, chunk_output_files, input_chunks, service)
        t = threading.Thread(
            target=run_chunk_wrapper,
            args=[chunks_in_flight, chunk_output_files, n, mutex, run_chunk,
                  [part_suffix, remote_home_dir, remote_index_dir, remote_work_dir, remote_username,
                   chunk_input_files, key_path, service, lazy_run]])
        t.start()
        chunk_threads.append(t)
    for ct in chunk_threads:
        ct.join()
        check_for_errors(mutex, chunk_output_files, input_chunks, service)
    assert None not in chunk_output_files
    concatenate_files(chunk_output_files, result_dir(output_file))
    async_handler.launch_aws_upload(RESULT_DIR + "/" + output_file, SAMPLE_S3_OUTPUT_PATH + "/")


def fetch_deuterostome_file(lock=threading.RLock()):  #pylint: disable=dangerous-default-value
    with lock:
        deuterostome_file_basename = os.path.basename(DEUTEROSTOME_TAXIDS)
        deuterostome_file = os.path.join(REF_DIR, deuterostome_file_basename)
        if not os.path.isfile(deuterostome_file):
            execute_command("aws s3 cp --quiet %s %s/" % (DEUTEROSTOME_TAXIDS, REF_DIR))
            write_to_log("downloaded deuterostome list")
        return deuterostome_file


def run_chunk_wrapper(chunks_in_flight, chunk_output_files, n, mutex, target, args):
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

def check_for_errors(mutex, chunk_output_files, input_chunks, what):
    with mutex:
        if "error" in chunk_output_files:
            # we alreay have per-chunk retries to address transient (non-deterministic) system issues;
            # if a chunk fails on all retries, it must be due to a deterministic & reproducible problem
            # with the chunk data or command options;  we should not even attempt the other chunks
            ei = chunk_output_files.index("error")
            raise Exception("All retries failed for {} chunk {}.".format(what, input_chunks[ei]))


def run_generate_unidentified_fasta(input_fa, output_fa):
    #TODO  remove annotated fasta intermediate file and replace > with : below
    subprocess.check_output("grep -A 1 '>NR::NT::' %s | sed '/^--$/d' > %s" % (input_fa, output_fa), shell=True)
    write_to_log("finished job")
    async_handler.launch_aws_upload(output_fa, SAMPLE_S3_OUTPUT_PATH + "/")


def get_alignment_version_s3_path():
    return os.path.join(SAMPLE_S3_OUTPUT_PATH, VERSION_OUT)


def fetch_input_and_replace_whitespace(input_filename, result):
    s3_input_path = os.path.join(SAMPLE_S3_INPUT_PATH, input_filename)
    s3_output_path = os.path.join(SAMPLE_S3_OUTPUT_PATH, input_filename) # usually the same as the input path
    cleaned_input_path = result_dir("nospace.%s" % input_filename)
    if not check_s3_file_presence(s3_input_path):
        result[0] = None
    else:
        execute_command("aws s3 cp --quiet {s3_input_path} - | sed 's/[[:blank:]]/{replacement}/g' > {cleaned_input_path}".format(
            replacement=";",
            s3_input_path=s3_input_path,
            cleaned_input_path=cleaned_input_path))
        result[0] = cleaned_input_path
        # This is extremely rare (debug/development only?) and we don't care if it succeeds
        if s3_input_path != s3_output_path:
            async_handler.launch_aws_upload(s3_input_path, s3_output_path)

def get_input_file_list():
    # Check existence of gsnap filter output
    if check_s3_file_presence(os.path.join(SAMPLE_S3_INPUT_PATH,
                                           EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1)):
        return [EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT1,
                EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT2,
                EXTRACT_UNMAPPED_FROM_GSNAP_SAM_OUT3]
    return [EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT1,
            EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT2,
            EXTRACT_UNMAPPED_FROM_BOWTIE_SAM_OUT3]

def fetch_and_clean_inputs():
    input_file_list = get_input_file_list()
    # Fetch inputs and remove tabs in parallel.
    cleaned_inputs = [
        ["ERROR"],
        ["ERROR"],
        ["ERROR"]
    ]
    input_fetcher_threads = [
        threading.Thread(target=fetch_input_and_replace_whitespace, args=[input_file_list[0], cleaned_inputs[0]]),
        threading.Thread(target=fetch_input_and_replace_whitespace, args=[input_file_list[1], cleaned_inputs[1]]),
        threading.Thread(target=fetch_input_and_replace_whitespace, args=[input_file_list[2], cleaned_inputs[2]])
    ]
    for ift in input_fetcher_threads:
        ift.start()
    for ift in input_fetcher_threads:
        ift.join()

    # flatten
    cleaned_inputs = [ci[0] for ci in cleaned_inputs]

    for i, ci in enumerate(cleaned_inputs):
        assert ci != "ERROR", "Error fetching input {}".format(i)
        assert ci == None or os.path.isfile(ci), "Local file {} not found after download of input {}".format(ci, i)

    assert cleaned_inputs[0] != None, "Input {} not found.  This input is mandatory.".format(input_file_list[0])

    assert (cleaned_inputs[1] == None) == (cleaned_inputs[2] == None), "Input {} is required when {} is given, and vice versa".format(input_file_list[1], input_file_list[2])

    cleaned_inputs = filter(None, cleaned_inputs)
    assert len(cleaned_inputs) in (1, 3)

    return cleaned_inputs


def run_stage2(lazy_run=True):
    # make local directories
    execute_command("mkdir -p %s %s %s %s %s" % (SAMPLE_DIR, FASTQ_DIR, RESULT_DIR, CHUNKS_RESULT_DIR, REF_DIR))

    # configure logger
    log_file = "%s/%s.%s.txt" % (RESULT_DIR, LOGS_OUT_BASENAME, AWS_BATCH_JOB_ID)
    configure_logger(log_file)

    # Initiate fetch of references
    lineage_paths = [None]
    accession2taxid_paths = [None]
    def fetch_references():
        lineage_paths[0] = fetch_reference(LINEAGE_SHELF)
        accession2taxid_paths[0] = fetch_reference(ACCESSION2TAXID)
    reference_fetcher_thread = threading.Thread(target=fetch_references)
    reference_fetcher_thread.start()

    cleaned_inputs = fetch_and_clean_inputs()

    gsnapl_input_files = [os.path.basename(f) for f in cleaned_inputs[:2]]
    merged_fasta = cleaned_inputs[-1]
    before_file_name_for_log = cleaned_inputs[0]
    before_file_type_for_log = "fasta_paired" if len(gsnapl_input_files) == 2 else "fasta" # unpaired

    # Import existing job stats
    stats = StatsFile(STATS_OUT, RESULT_DIR, SAMPLE_S3_INPUT_PATH, SAMPLE_S3_OUTPUT_PATH)
    stats.load_from_s3()

    thread_success = {}
    thread_success_lock = threading.RLock()

    # subsample if specified
    if SUBSAMPLE:
        target_n_reads = int(SUBSAMPLE)
        if len(gsnapl_input_files) == 2:
            subsampled_gsnapl_input_files, subsampled_merged_fasta = subsample_paired_fastas(gsnapl_input_files, os.path.basename(merged_fasta), target_n_reads)
        else:
            subsampled_merged_fasta = subsample_single_fasta(gsnapl_input_files[0], target_n_reads)
            subsampled_gsnapl_input_files = [subsampled_merged_fasta]
        gsnapl_input_files = subsampled_gsnapl_input_files
        merged_fasta = result_dir(subsampled_merged_fasta)
        for i, f in enumerate(gsnapl_input_files):
            async_handler.launch_aws_upload(result_dir(f), SAMPLE_S3_OUTPUT_PATH + "/")
        if len(gsnapl_input_files) == 2:
            async_handler.launch_aws_upload(merged_fasta, SAMPLE_S3_OUTPUT_PATH + "/")

    deuterostome_fetcher = threading.Thread(target=fetch_deuterostome_file)
    deuterostome_fetcher.start()


    def run_gsnap():
        # run gsnap remotely
        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "GSNAPL",
             "version_file_s3": GSNAP_VERSION_FILE_S3,
             "output_version_file": result_dir(VERSION_OUT)})
        run_and_log_s3(
            logparams,
            [result_dir(MULTIHIT_GSNAPL_OUT)],
            lazy_run,
            SAMPLE_S3_OUTPUT_PATH,
            run_remotely,
            gsnapl_input_files,
            "gsnap",
            lazy_run)


        reference_fetcher_thread.join()
        call_hits_m8(
            result_dir(MULTIHIT_GSNAPL_OUT),
            lineage_paths[0],
            accession2taxid_paths[0],
            result_dir(DEDUP_MULTIHIT_GSNAPL_OUT),
            result_dir(SUMMARY_MULTIHIT_GSNAPL_OUT)
        )
        stats.count_reads("run_gsnapl_remotely",
                          before_filename=before_file_name_for_log,
                          before_filetype=before_file_type_for_log,
                          after_filename=result_dir(DEDUP_MULTIHIT_GSNAPL_OUT),
                          after_filetype="m8")

        deuterostome_path = None
        if not SKIP_DEUTERO_FILTER:
            deuterostome_fetcher.join()
            deuterostome_path = fetch_deuterostome_file()

        generate_taxon_count_json_from_m8(result_dir(DEDUP_MULTIHIT_GSNAPL_OUT),
                                          result_dir(SUMMARY_MULTIHIT_GSNAPL_OUT),
                                          'raw',
                                          'NT',
                                          lineage_paths[0],
                                          deuterostome_path,
                                          stats.get_total_reads(),
                                          stats.get_remaining_reads(),
                                          result_dir(MULTIHIT_NT_JSON_OUT))

        with thread_success_lock:
            thread_success["gsnap"] = True


    def run_rapsearch2():
        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "RAPSearch2",
             "version_file_s3": RAPSEARCH_VERSION_FILE_S3, "output_version_file": result_dir(VERSION_OUT)})
        run_and_log_s3(
            logparams,
            [result_dir(MULTIHIT_RAPSEARCH_OUT)],
            lazy_run,
            SAMPLE_S3_OUTPUT_PATH,
            run_remotely,
            [merged_fasta],
            "rapsearch",
            lazy_run)

        reference_fetcher_thread.join()
        call_hits_m8(
            result_dir(MULTIHIT_RAPSEARCH_OUT),
            lineage_paths[0],
            accession2taxid_paths[0],
            result_dir(DEDUP_MULTIHIT_RAPSEARCH_OUT),
            result_dir(SUMMARY_MULTIHIT_RAPSEARCH_OUT)
        )
        stats.count_reads("run_rapsearch2_remotely",
                          before_filename=result_dir(merged_fasta),
                          before_filetype="fasta",
                          after_filename=result_dir(DEDUP_MULTIHIT_RAPSEARCH_OUT),
                          after_filetype="m8")

        deuterostome_path = None
        if not SKIP_DEUTERO_FILTER:
            deuterostome_fetcher.join()
            deuterostome_path = fetch_deuterostome_file()

        # PRODUCE NEW MULTIHIT NR OUTPUT
        generate_taxon_count_json_from_m8(result_dir(DEDUP_MULTIHIT_RAPSEARCH_OUT),
                                          result_dir(SUMMARY_MULTIHIT_RAPSEARCH_OUT),
                                          'log10',
                                          'NR',
                                          lineage_paths[0],
                                          deuterostome_path,
                                          stats.get_total_reads(),
                                          stats.get_remaining_reads(),
                                          result_dir(MULTIHIT_NR_JSON_OUT))

        with thread_success_lock:
            thread_success["rapsearch2"] = True


    def run_additional_steps():

        # COMBINE NEW MULTIHIT NT AND NR OUTPUTS
        combine_pipeline_output_json(result_dir(MULTIHIT_NT_JSON_OUT),
                                     result_dir(MULTIHIT_NR_JSON_OUT),
                                     result_dir(MULTIHIT_COMBINED_JSON_OUT),
                                     stats)
        # Keep this warning there for a month or so, long enough for obselete webapps to retire.
        # We have to do this because an obsolete webapp will keep the gsnap tier from scaling in if it doesn't see this file.
        execute_command("echo This file is deprecated since pipeline version 1.5.  Please use {new} instead. > {deprecated}".format(
            deprecated=result_dir(DEPRECATED_BOOBYTRAPPED_COMBINED_JSON_OUT),
            new=MULTIHIT_COMBINED_JSON_OUT))
        src = RESULT_DIR + "/" + DEPRECATED_BOOBYTRAPPED_COMBINED_JSON_OUT
        dst = SAMPLE_S3_OUTPUT_PATH + "/"
        async_handler.launch_aws_upload(src, dst)

        with thread_success_lock:
            thread_success["additional_steps"] = True


    def run_annotation_steps():
        # Annotate fasta with accessions
        annotate_fasta_with_accessions(merged_fasta,
                                       result_dir(DEDUP_MULTIHIT_GSNAPL_OUT),
                                       result_dir(DEDUP_MULTIHIT_RAPSEARCH_OUT),
                                       result_dir(ACCESSION_ANNOTATED_FASTA))

        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "generate FASTA of unidentified reads"})
        run_and_log_eager(
            logparams,
            [result_dir(UNIDENTIFIED_FASTA_OUT)],
            run_generate_unidentified_fasta,
            result_dir(ACCESSION_ANNOTATED_FASTA),
            result_dir(UNIDENTIFIED_FASTA_OUT))
        stats.count_reads("run_generate_unidentified_fasta",
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

    t_gsnap.start()
    t_rapsearch2.start()
    t_gsnap.join()
    assert thread_succeeded("gsnap")
    t_rapsearch2.join()
    assert thread_succeeded("rapsearch2")

    t_additional_steps.start()
    t_annotation.start()
    t_additional_steps.join()
    assert thread_succeeded("additional_steps")
    t_annotation.join()
    assert thread_succeeded("annotation")

    async_handler.wait_on_all()

    # copy log file -- after work is done
    write_to_log("Non-host alignment complete")
    stats.save_to_s3()
    upload_log_file(SAMPLE_S3_OUTPUT_PATH)
