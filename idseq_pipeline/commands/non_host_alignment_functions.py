import os
import subprocess
import json
import shelve
import re
import time
import math
import threading
import shutil
import traceback
import random
import multiprocessing
from .common import * #pylint: disable=wildcard-import


# that's per job;  by default in early 2018 jobs are subsampled to <= 100 chunks
MAX_CHUNKS_IN_FLIGHT_GSNAP = 64
MAX_CHUNKS_IN_FLIGHT_RAPSEARCH = 64
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

# output files
EXTRACT_UNMAPPED_FROM_SAM_OUT1 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.1.fasta'
EXTRACT_UNMAPPED_FROM_SAM_OUT2 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.2.fasta'
EXTRACT_UNMAPPED_FROM_SAM_OUT3 = 'unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.merged.fasta'
GSNAPL_OUT = 'gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
GSNAPL_DEDUP_OUT = 'dedup.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
ANNOTATE_GSNAPL_M8_WITH_TAXIDS_OUT = 'taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
GENERATE_TAXID_ANNOTATED_FASTA_FROM_M8_OUT = 'taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.fasta'
RAPSEARCH2_OUT = 'rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
FILTER_DEUTEROSTOMES_FROM_NT_M8_OUT = 'filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
NT_M8_TO_TAXID_COUNTS_FILE_OUT = 'counts.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.csv'
NT_TAXID_COUNTS_TO_JSON_OUT = 'counts.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.json'
ANNOTATE_RAPSEARCH2_M8_WITH_TAXIDS_OUT = 'taxids.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
GENERATE_TAXID_ANNOTATED_FASTA_FROM_RAPSEARCH2_M8_OUT = 'taxids.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.fasta'
FILTER_DEUTEROSTOMES_FROM_NR_M8_OUT = 'filter.deuterostomes.taxids.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
NR_M8_TO_TAXID_COUNTS_FILE_OUT = 'counts.filter.deuterostomes.taxids.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.csv'
NR_TAXID_COUNTS_TO_JSON_OUT = 'counts.filter.deuterostomes.taxids.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.json'
UNIDENTIFIED_FASTA_OUT = 'unidentified.fasta'
COMBINED_JSON_OUT = 'idseq_web_sample.json'
LOGS_OUT_BASENAME = 'log'
STATS_OUT = 'stats.json'
VERSION_OUT = 'versions.json'

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
CHUNKS_RESULT_DIR = os.path.join(RESULT_DIR, "chunks")
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

# target outputs by task
TARGET_OUTPUTS = {"run_gsnapl_remotely": [os.path.join(RESULT_DIR, GSNAPL_DEDUP_OUT)],
                  "run_annotate_m8_with_taxids__1": [os.path.join(RESULT_DIR, ANNOTATE_GSNAPL_M8_WITH_TAXIDS_OUT)],
                  "run_generate_taxid_annotated_fasta_from_m8__1": [os.path.join(RESULT_DIR, GENERATE_TAXID_ANNOTATED_FASTA_FROM_M8_OUT)],
                  "run_filter_deuterostomes_from_m8__1": [os.path.join(RESULT_DIR, FILTER_DEUTEROSTOMES_FROM_NT_M8_OUT)],
                  "run_generate_taxid_outputs_from_m8__1": [os.path.join(RESULT_DIR, NT_TAXID_COUNTS_TO_JSON_OUT)],
                  "run_rapsearch2_remotely": [os.path.join(RESULT_DIR, RAPSEARCH2_OUT)],
                  "run_annotate_m8_with_taxids__2": [os.path.join(RESULT_DIR, ANNOTATE_RAPSEARCH2_M8_WITH_TAXIDS_OUT)],
                  "run_generate_taxid_annotated_fasta_from_m8__2": [os.path.join(RESULT_DIR, GENERATE_TAXID_ANNOTATED_FASTA_FROM_RAPSEARCH2_M8_OUT)],
                  "run_filter_deuterostomes_from_m8__2": [os.path.join(RESULT_DIR, FILTER_DEUTEROSTOMES_FROM_NR_M8_OUT)],
                  "run_generate_taxid_outputs_from_m8__2": [os.path.join(RESULT_DIR, NR_TAXID_COUNTS_TO_JSON_OUT)],
                  "run_combine_json_outputs": [os.path.join(RESULT_DIR, COMBINED_JSON_OUT)],
                  "run_generate_unidentified_fasta": [os.path.join(RESULT_DIR, UNIDENTIFIED_FASTA_OUT)]}

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

# convenience functions
def count_lines_in_paired_files(input_files):
    known_nlines = None
    for input_file in input_files:
        nlines = int(execute_command_with_output("wc -l %s" % input_file).strip().split()[0])
        if known_nlines != None:
            assert nlines == known_nlines, "Mismatched line counts in supposedly paired files: {}".format(input_files)
        known_nlines = nlines
    return nlines

def subsample_single_fasta(input_file, records_to_keep, type_, output_file):
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

def subsample_fastas(input_files_basenames, merged_file_basename, target_n_reads):
    input_files = [os.path.join(RESULT_DIR, f) for f in input_files_basenames]
    merged_file = os.path.join(RESULT_DIR, merged_file_basename)
    total_records = count_lines_in_paired_files(input_files) // 2 # each fasta record spans 2 lines
    write_to_log("total read pairs: %d" % total_records)
    write_to_log("target read pairs: %d" % target_n_reads)
    # note: target_n_reads and total_records really refer to numbers of read PAIRS
    if total_records <= target_n_reads:
        return input_files_basenames, merged_file_basename
    subsample_prefix = "subsample_%d" % target_n_reads
    records_to_keep = set(random.sample(xrange(total_records), target_n_reads))
    subsampled_files = []
    known_kept_read_ids = set()
    kept_reads_count = 0
    # subsample the paired files and record read IDs kept
    for input_file in input_files:
        input_dir = os.path.split(input_file)[0]
        input_basename = os.path.split(input_file)[1]
        output_basename = "%s.%s" % (subsample_prefix, input_basename)
        output_file = os.path.join(input_dir, output_basename)
        kept_read_ids, kept_reads_count = subsample_single_fasta(input_file, records_to_keep, "record_indices", output_file)
        assert kept_reads_count == target_n_reads
        subsampled_files.append(output_basename)
        known_kept_read_ids |= kept_read_ids
    # subsample the merged file to the same read IDs
    input_dir = os.path.split(merged_file)[0]
    input_basename = os.path.split(merged_file)[1]
    subsampled_merged_file = "%s.%s" % (subsample_prefix, input_basename)
    output_file = os.path.join(input_dir, subsampled_merged_file)
    _, kept_reads_count = subsample_single_fasta(merged_file, known_kept_read_ids, "read_ids", output_file)
    num_desired = target_n_reads * len(input_files_basenames)
    if kept_reads_count != num_desired:
        # TODO Need a way to surface these warnings somehow in the webapp.
        write_to_log("ERROR:  Improperly paired reads.  Total reads in sampled merged fasta: {krc}".format(krc=kept_reads_count))
    assert kept_reads_count == num_desired
    return subsampled_files, subsampled_merged_file


def concatenate_files(file_list, output_file):
    with open(output_file, 'wb') as outf:
        for f in file_list:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, outf)

def remove_annotation(read_id):
    result = re.sub(r'NT:(.*?):', '', read_id)
    result = re.sub(r'NR:(.*?):', '', result)
    return result

def generate_taxid_annotated_fasta_from_m8(input_fasta_file, m8_file, output_fasta_file, annotation_prefix):
    '''Tag reads based on the m8 output'''
    # Example:  generate_annotated_fasta_from_m8('filter.unmapped.merged.fasta',
    #  'bowtie.unmapped.star.gsnapl-nt-k16.m8', 'NT-filter.unmapped.merged.fasta', 'NT')
    # Construct the m8_hash
    read_to_accession_id = {}
    with open(m8_file, 'rb') as m8f:
        for line in m8f:
            if line[0] == '#':
                continue
            parts = line.split("\t")
            log_corrupt(len(parts) < 12, m8_file, line)
            read_name = parts[0]
            read_name_parts = read_name.split("/")
            if len(read_name_parts) > 1:
                output_read_name = read_name_parts[0] + '/' + read_name_parts[-1]
            else:
                output_read_name = read_name
            accession_id = parts[1]
            read_to_accession_id[output_read_name] = accession_id
    # Go through the input_fasta_file to get the results and tag reads
    input_fasta_f = open(input_fasta_file, 'rb')
    output_fasta_f = open(output_fasta_file, 'wb')
    sequence_name = input_fasta_f.readline()
    sequence_data = input_fasta_f.readline()
    while sequence_name and sequence_data:
        read_id = sequence_name.rstrip().lstrip('>')
        accession = read_to_accession_id.get(remove_annotation(read_id), '')
        new_read_name = annotation_prefix + ':' + accession + ':' + read_id
        output_fasta_f.write(">%s\n" % new_read_name)
        output_fasta_f.write(sequence_data)
        sequence_name = input_fasta_f.readline()
        sequence_data = input_fasta_f.readline()
    input_fasta_f.close()
    output_fasta_f.close()


def deduplicate_m8(input_m8, output_m8):
    outf = open(output_m8, "wb")
    previous_read_name = ''
    with open(input_m8, "rb") as m8f:
        for line in m8f:
            if line[0] == '#':
                continue
            parts = line.split("\t")
            read_name = parts[0] # Example: HWI-ST640:828:H917FADXX:2:1108:8883:88679/1/1'
            if read_name == previous_read_name:
                continue
            outf.write(line)
            previous_read_name = read_name
    outf.close()


def log_corrupt(is_corrupt, m8_file, line):
    if is_corrupt:
        write_to_log(m8_file + " is corrupt at line:\n" + line + "\n----> delete it and its corrupt ancestors before restarting run")
        raise AssertionError


def generate_tax_counts_from_m8(m8_file, e_value_type, output_file, lineage_map):
    taxid_count_map = {}
    taxid_percent_identity_map = {}
    taxid_alignment_length_map = {}
    taxid_e_value_map = {}
    read_to_taxid = {}

    with open(m8_file, 'rb') as m8f:
        for line in m8f:
            # Get taxid:
            line_columns = line.split("\t")
            # If file is corrupt, log it:
            log_corrupt(len(line_columns) < 12, m8_file, line)
            # Otherwise continue processing:
            read_id_column = line_columns[0]
            taxid = (read_id_column.split("taxid"))[1].split(":")[0]
            species_taxid, genus_taxid, family_taxid = lineage_map.get(taxid, ("-100", "-200", "-300"))

            raw_read_id = extract_m8_readid(read_id_column)

            # Get alignment quality metrics from m8 format (blast format 8):
            #   %id, alignment length, mismatches, gap openings, query start, query end,
            #   subject start, subject end, E value (log10 if rapsearch2 output), bit score
            # E value is a negative power of 10. GSNAPL outputs raw e-value, RAPSearch2 outputs log10(e-value).
            # Whenever we use "e_value" it refers to log10(e-value), which is easier to handle.
            percent_identity = float(line_columns[2])
            alignment_length = float(line_columns[3])
            e_value = float(line_columns[10])
            if e_value_type != 'log10':
                e_value = math.log10(e_value)
            taxid_count_map[taxid] = taxid_count_map.get(taxid, 0) + 1
            taxid_percent_identity_map[taxid] = taxid_percent_identity_map.get(taxid, 0) + percent_identity
            taxid_alignment_length_map[taxid] = taxid_alignment_length_map.get(taxid, 0) + alignment_length
            taxid_e_value_map[taxid] = taxid_e_value_map.get(taxid, 0) + e_value

            # Keep mapping of read ids to taxids in order to determine pair concordance later:
            read_to_taxid[raw_read_id] = (species_taxid, genus_taxid, family_taxid)

    # Write results:
    with open(output_file, 'w') as f:
        for taxid, count in taxid_count_map.iteritems():
            avg_percent_identity = taxid_percent_identity_map[taxid] / count
            avg_alignment_length = taxid_alignment_length_map[taxid] / count
            avg_e_value = taxid_e_value_map[taxid] / count
            f.write(",".join([str(taxid), str(count), str(avg_percent_identity), str(avg_alignment_length), str(avg_e_value) + '\n']))
    # Determine pair concordance:
    return check_pair_concordance(read_to_taxid)


def check_pair_concordance(read_to_taxid):
    species_taxid_concordance_map = {}
    genus_taxid_concordance_map = {}
    family_taxid_concordance_map = {}
    read_1_ids = [read_id for read_id in read_to_taxid.keys() if "/1" in read_id]
    for read_1_id in read_1_ids:
        read_2_id = read_1_id.replace("/1", "/2", 1)
        if read_2_id not in read_to_taxid:
            continue
        species_taxid_1, genus_taxid_1, family_taxid_1 = read_to_taxid[read_1_id]
        species_taxid_2, genus_taxid_2, family_taxid_2 = read_to_taxid[read_2_id]
        if species_taxid_1 == species_taxid_2: # add both reads to the concordance count
            species_taxid_concordance_map[species_taxid_1] = species_taxid_concordance_map.get(species_taxid_1, 0) + 2
        if genus_taxid_1 == genus_taxid_2:
            genus_taxid_concordance_map[genus_taxid_1] = genus_taxid_concordance_map.get(genus_taxid_1, 0) + 2
        if family_taxid_1 == family_taxid_2:
            family_taxid_concordance_map[family_taxid_1] = family_taxid_concordance_map.get(family_taxid_1, 0) + 2
    return species_taxid_concordance_map, genus_taxid_concordance_map, family_taxid_concordance_map

def generate_json_from_taxid_counts(taxidCountsInputPath, jsonOutputPath, countType, lineage_map,
                                    species_total_concordant, genus_total_concordant, family_total_concordant, stats):
    total_reads = stats.get_total_reads()
    taxon_counts_attributes = []
    remaining_reads = stats.get_remaining_reads()

    species_to_count = {}
    species_to_percent_identity = {}
    species_to_alignment_length = {}
    species_to_e_value = {}
    species_to_species_level_concordance = {}
    species_to_genus_level_concordance = {}
    species_to_family_level_concordance = {}
    with open(taxidCountsInputPath) as f:
        for line in f:
            tok = line.rstrip().split(",")
            taxid = tok[0]
            count = float(tok[1])
            percent_identity = float(tok[2])
            alignment_length = float(tok[3])
            e_value = float(tok[4])
            species_taxid, genus_taxid, family_taxid = lineage_map.get(taxid, ("-100", "-200", "-300"))
            species_to_count[species_taxid] = species_to_count.get(species_taxid, 0) + count
            species_to_percent_identity[species_taxid] = species_to_percent_identity.get(species_taxid, 0) + count * percent_identity
            species_to_alignment_length[species_taxid] = species_to_alignment_length.get(species_taxid, 0) + count * alignment_length
            species_to_e_value[species_taxid] = species_to_e_value.get(species_taxid, 0) + count * e_value
            species_to_species_level_concordance[species_taxid] = species_total_concordant.get(species_taxid, 0)
            species_to_genus_level_concordance[species_taxid] = genus_total_concordant.get(genus_taxid, 0)
            species_to_family_level_concordance[species_taxid] = family_total_concordant.get(family_taxid, 0)

    for taxid in species_to_count:
        count = species_to_count[taxid]
        avg_percent_identity = species_to_percent_identity[taxid] / count
        avg_alignment_length = species_to_alignment_length[taxid] / count
        avg_e_value = species_to_e_value[taxid] / count
        taxon_counts_attributes.append({"tax_id": taxid,
                                        "tax_level": TAX_LEVEL_SPECIES,
                                        "count": count,
                                        "percent_identity": avg_percent_identity,
                                        "alignment_length": avg_alignment_length,
                                        "e_value": avg_e_value,
                                        "count_type": countType,
                                        "percent_concordant": (100.0 * species_to_species_level_concordance[taxid]) / count,
                                        # Not very elegant, but until such time as we propagate alignment information at the level of
                                        # individual reads to the web app's database, we have to do the concordance aggregation here:
                                        "species_total_concordant": species_to_species_level_concordance[taxid],
                                        "genus_total_concordant": species_to_genus_level_concordance[taxid],
                                        "family_total_concordant": species_to_family_level_concordance[taxid]})

    output_dict = {
        "pipeline_output": {
            "total_reads": total_reads,
            "remaining_reads": remaining_reads,
            "taxon_counts_attributes": taxon_counts_attributes
        }
    }
    with open(jsonOutputPath, 'wb') as outf:
        json.dump(output_dict, outf)


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


def generate_taxid_annotated_m8(input_m8, output_m8, accession2taxid_db):
    accession2taxid_dict = shelve.open(accession2taxid_db)
    outf = open(output_m8, "wb")
    with open(input_m8, "rb") as m8f:
        for line in m8f:
            if line[0] == '#':
                continue
            parts = line.split("\t")
            log_corrupt(len(parts) < 12, input_m8, line)
            _read_name = parts[0] # Example: HWI-ST640:828:H917FADXX:2:1108:8883:88679/1/1',
            accession_id = parts[1] # Example: CP000671.1',
            accession_id_short = accession_id.split(".")[0]
            new_line = "taxid" + accession2taxid_dict.get(accession_id_short, "NA") + ":" + line
            outf.write(new_line)
    outf.close()


def generate_merged_fasta(input_files, output_file):
    with open(output_file, 'w') as outfile:
        for fname in input_files:
            idx = input_files.index(fname) + 1
            with open(fname) as infile:
                for line in infile:
                    if line.startswith(">") and not "/" in line:
                        suffix = "/" + str(idx)
                    else:
                        suffix = ""
                    outfile.write(line.rstrip() + suffix + "\n")


def read_file_into_set(file_name):
    with open(file_name) as f:
        S = set(x.rstrip() for x in f)
    S.discard('')
    return S


def filter_deuterostomes_from_m8(input_m8, output_m8, deuterostome_file, cache={}, lock=threading.RLock()):  #pylint: disable=dangerous-default-value
    with lock:
        if deuterostome_file not in cache:
            cache[deuterostome_file] = read_file_into_set(deuterostome_file)
        taxids_toremove = cache[deuterostome_file]
    output_f = open(output_m8, 'wb')
    with open(input_m8, "rb") as input_f:
        for line in input_f:
            # Every line should start with "taxid<number>:...", for example
            #    "taxid9606:NB501961:14:HM7TLBGX2:1:12104:15431:..."
            # Exceptional case:  If a line doesn't match that pattern, we output the line anyway.
            taxid = None
            if line.startswith("taxid"):
                taxid = line.split(":", 1)[0][5:]
            if taxid not in taxids_toremove:
                output_f.write(line)
    output_f.close()


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
        input_file_full_local_path = os.path.join(RESULT_DIR, input_file)
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


def remove_whitespace_from_files(input_files, replacement, output_files):
    for idx, input_file in enumerate(input_files):
        output_file = output_files[idx]
        execute_command("sed 's/[[:blank:]]/%s/g' %s > %s" % (replacement, input_file, output_file))


def clean_direct_gsnapl_input(fastq_files):
    # unzip files if necessary
    if ".gz" in FILE_TYPE:
        subprocess.check_output(" ".join(["gunzip", "-f"] + fastq_files), shell=True)
        unzipped_files = [os.path.splitext(f)[0] for f in fastq_files]
    else:
        unzipped_files = fastq_files
    # convert to fasta if necessary
    file_type_trimmed = FILE_TYPE.split(".gz")[0]
    if file_type_trimmed == "fastq":
        file_prefixes = [os.path.splitext(f)[0] for f in unzipped_files]
        for file_prefix in file_prefixes:
            execute_command("sed -n '1~4s/^@/>/p;2~4p' <%s.fastq >%s.fasta" % (file_prefix, file_prefix))
        cleaned_files = [f + ".fasta" for f in file_prefixes]
    else:
        cleaned_files = unzipped_files
    # generate file type for log
    file_type_for_log = "fasta"
    if len(fastq_files) == 2:
        file_type_for_log += "_paired"
    # copy files to S3
    for f in cleaned_files:
        execute_command("aws s3 cp --quiet %s %s/" % (f, SAMPLE_S3_OUTPUT_PATH))
    return cleaned_files, file_type_for_log


def interpret_min_column_number_string(min_column_number_string, correct_number_of_output_columns, try_number):
    if min_column_number_string:
        min_column_number = float(min_column_number_string)
        write_to_log("Try no. %d: Smallest number of columns observed in any line was %d" % (try_number, min_column_number))
    else:
        write_to_log("Try no. %d: No hits" % try_number)
        min_column_number = correct_number_of_output_columns
    return min_column_number


def run_gsnapl_chunk(part_suffix, remote_home_dir, remote_index_dir, remote_work_dir, remote_username,
                     input_files, key_path, lazy_run):
    chunk_id = input_files[0].split(part_suffix)[-1]
    outfile_basename = 'gsnapl-out' + part_suffix + chunk_id
    dedup_outfile_basename = 'dedup-' + outfile_basename
    local_outfile = CHUNKS_RESULT_DIR + "/" + outfile_basename
    remote_outfile = os.path.join(remote_work_dir, outfile_basename)
    commands = "mkdir -p %s;" % remote_work_dir
    for input_fa in input_files:
        commands += "aws s3 cp --quiet %s/%s %s/ ; " % \
                 (SAMPLE_S3_OUTPUT_CHUNKS_PATH, input_fa, remote_work_dir)
    commands += " ".join([remote_home_dir+'/bin/gsnapl',
                          '-A', 'm8', '--batch=0', '--use-shared-memory=0',
                          '--gmap-mode=none', '--npaths=1', '--ordered',
                          '-t', '32',
                          '--maxsearch=5', '--max-mismatches=40',
                          '-D', remote_index_dir, '-d', 'nt_k16']
                         + [remote_work_dir+'/'+input_fa for input_fa in input_files]
                         + ['> '+remote_outfile, ';'])
    # Also produce gsnap output with mutliple hits per read
    multihit_local_outfile = CHUNKS_RESULT_DIR + "/multihit-" + outfile_basename
    multihit_remote_outfile = 'multihit-' + remote_outfile
    commands += " ".join([remote_home_dir+'/bin/gsnapl',
                          '-A', 'm8', '--batch=0', '--use-shared-memory=0',
                          '--gmap-mode=none', '--npaths=1000', '--ordered',
                          '-t', '32',
                          '--maxsearch=1000', '--max-mismatches=40',
                          '-D', remote_index_dir, '-d', 'nt_k16']
                         + [remote_work_dir+'/'+input_fa for input_fa in input_files]
                         + ['> ' + multihit_remote_outfile, ';'])

    if not lazy_run or not fetch_lazy_result(os.path.join(SAMPLE_S3_OUTPUT_CHUNKS_PATH, dedup_outfile_basename), CHUNKS_RESULT_DIR):
        correct_number_of_output_columns = 12
        min_column_number = 0
        max_tries = 2
        try_number = 1
        while min_column_number != correct_number_of_output_columns and try_number <= max_tries:
            write_to_log("waiting for gsnap server for chunk {}".format(chunk_id))
            gsnapl_instance_ip = wait_for_server_ip('gsnap', key_path, remote_username, ENVIRONMENT, GSNAPL_MAX_CONCURRENT, chunk_id)
            write_to_log("starting alignment for chunk %s on machine %s" % (chunk_id, gsnapl_instance_ip))
            execute_command(remote_command(commands, key_path, remote_username, gsnapl_instance_ip))
            # check if every row has correct number of columns (12) in the output file on the remote machine
            verification_command = "awk '{print NF}' %s | sort -nu | head -n 1" % remote_outfile
            min_column_number_string = execute_command_with_output(remote_command(verification_command, key_path, remote_username, gsnapl_instance_ip))
            min_column_number = interpret_min_column_number_string(min_column_number_string, correct_number_of_output_columns, try_number)
            try_number += 1
        # move output from remote machine to s3
        assert min_column_number == correct_number_of_output_columns, "Chunk %s output corrupt; not copying to S3. Re-start pipeline to try again." % chunk_id
        with iostream:
            execute_command(scp(key_path, remote_username, gsnapl_instance_ip, remote_outfile, local_outfile))
            execute_command(scp(key_path, remote_username, gsnapl_instance_ip, multihit_remote_outfile, multihit_local_outfile))
        # Deduplicate m8 input. Sometimes GSNAPL outputs multiple consecutive lines for same original read and same accession id. Count functions expect only 1 (top hit).
        execute_command("echo '' >> %s;" % local_outfile) # add a blank line at the end of the file so S3 copy doesn't fail if output is empty
        deduplicate_m8(os.path.join(CHUNKS_RESULT_DIR, outfile_basename), os.path.join(CHUNKS_RESULT_DIR, dedup_outfile_basename))
        with iostream:
            execute_command("aws s3 cp --quiet %s/%s %s/" % (CHUNKS_RESULT_DIR, dedup_outfile_basename, SAMPLE_S3_OUTPUT_CHUNKS_PATH))
    with iostream:
        execute_command("sed -i '$ {/^$/d;}' %s" % os.path.join(CHUNKS_RESULT_DIR, dedup_outfile_basename)) # remove blank line from end of file
    write_to_log("finished alignment for chunk %s" % chunk_id)
    return os.path.join(CHUNKS_RESULT_DIR, dedup_outfile_basename)



def run_gsnapl_remotely(input_files, lazy_run):
    output_file = os.path.join(SAMPLE_S3_OUTPUT_PATH, GSNAPL_DEDUP_OUT)
    key_path = fetch_key(ENVIRONMENT)
    remote_username = "ubuntu"
    remote_home_dir = "/home/%s" % remote_username
    remote_work_dir = "%s/batch-pipeline-workdir/%s" % (remote_home_dir, SAMPLE_NAME)
    remote_index_dir = "%s/share" % remote_home_dir
    # split file:
    chunk_nlines = 2*GSNAPL_CHUNK_SIZE
    part_suffix, input_chunks = chunk_input(input_files, chunk_nlines, GSNAPL_CHUNK_SIZE)
    # process chunks:
    chunk_output_files = [None] * len(input_chunks)
    chunk_threads = []
    mutex = threading.RLock()
    iii = list(enumerate(input_chunks))
    random.shuffle(iii)
    chunks_in_flight = chunks_in_flight_gsnap
    for n, chunk_input_files in iii:
        chunks_in_flight.acquire()
        check_for_errors(mutex, chunk_output_files, input_chunks, "gsnap")
        t = threading.Thread(
            target=run_chunk_wrapper,
            args=[chunks_in_flight, chunk_output_files, n, mutex, run_gsnapl_chunk,
                  [part_suffix, remote_home_dir, remote_index_dir, remote_work_dir, remote_username,
                   chunk_input_files, key_path, lazy_run]])
        t.start()
        chunk_threads.append(t)
    for ct in chunk_threads:
        ct.join()
        check_for_errors(mutex, chunk_output_files, input_chunks, "gsnap")
    assert None not in chunk_output_files
    # merge output chunks:
    with iostream:
        concatenate_files(chunk_output_files, os.path.join(RESULT_DIR, GSNAPL_DEDUP_OUT))
        execute_command("aws s3 cp --quiet %s/%s %s" % (RESULT_DIR, GSNAPL_DEDUP_OUT, output_file))
        # Also copy multihit outputs
        multihit_chunk_output_files = [f.replace("dedup", "multihit") for f in chunk_output_files]
        concatenate_files(multihit_chunk_output_files, os.path.join(RESULT_DIR, MULTIHIT_GSNAPL_OUT))
        execute_command("aws s3 cp --quiet %s/%s %s" % (RESULT_DIR, MULTIHIT_GSNAPL_OUT, os.path.join(SAMPLE_S3_OUTPUT_PATH, MULTIHIT_GSNAPL_OUT)))



def run_annotate_m8_with_taxids(input_m8, output_m8):
    accession2taxid_path = fetch_reference(ACCESSION2TAXID)
    p = multiprocessing.Process(target=generate_taxid_annotated_m8, args=[input_m8, output_m8, accession2taxid_path])
    p.start()
    p.join()
    if p.exitcode != 0:
        raise Exception("Failed generate_taxid_annotated_m8 on {}".format(input_m8))
    write_to_log("finished annotation")
    # move the output back to S3
    execute_command("aws s3 cp --quiet %s %s/" % (output_m8, SAMPLE_S3_OUTPUT_PATH))


def fetch_deuterostome_file(lock=threading.RLock()):  #pylint: disable=dangerous-default-value
    with lock:
        deuterostome_file_basename = os.path.basename(DEUTEROSTOME_TAXIDS)
        deuterostome_file = os.path.join(REF_DIR, deuterostome_file_basename)
        if not os.path.isfile(deuterostome_file):
            execute_command("aws s3 cp --quiet %s %s/" % (DEUTEROSTOME_TAXIDS, REF_DIR))
            write_to_log("downloaded deuterostome list")
        return deuterostome_file


def run_filter_deuterostomes_from_m8(input_m8, output_m8):
    deuterostome_file = fetch_deuterostome_file()
    filter_deuterostomes_from_m8(input_m8, output_m8, deuterostome_file)
    write_to_log("finished job")
    # move the output back to S3
    execute_command("aws s3 cp --quiet %s %s/" % (output_m8, SAMPLE_S3_OUTPUT_PATH))


def run_generate_taxid_annotated_fasta_from_m8(input_m8, input_fasta,
                                               output_fasta, annotation_prefix):
    generate_taxid_annotated_fasta_from_m8(input_fasta, input_m8, output_fasta, annotation_prefix)
    write_to_log("finished job")
    # move the output back to S3
    execute_command("aws s3 cp --quiet %s %s/" % (output_fasta, SAMPLE_S3_OUTPUT_PATH))




def run_rapsearch_chunk(part_suffix, _remote_home_dir, remote_index_dir, remote_work_dir, remote_username,
                        input_fasta, key_path, lazy_run):
    chunk_id = input_fasta.split(part_suffix)[-1]
    commands = "mkdir -p %s;" % remote_work_dir
    commands += "aws s3 cp --quiet %s/%s %s/ ; " % \
                 (SAMPLE_S3_OUTPUT_CHUNKS_PATH, input_fasta, remote_work_dir)
    input_path = remote_work_dir + '/' + input_fasta
    outfile_basename = 'rapsearch2-out' + part_suffix + chunk_id + '.m8'
    output_path = os.path.join(remote_work_dir, outfile_basename)
    commands += " ".join(['/usr/local/bin/rapsearch',
                          '-d', remote_index_dir+'/nr_rapsearch',
                          '-e', '-6',
                          '-l', '10',
                          '-a', 'T',
                          '-b', '0',
                          '-v', '1',
                          '-z', '24',
                          '-q', input_path,
                          '-o', output_path[:-3],
                          ';'])
    if not lazy_run or not fetch_lazy_result(os.path.join(SAMPLE_S3_OUTPUT_CHUNKS_PATH, outfile_basename), CHUNKS_RESULT_DIR):
        correct_number_of_output_columns = 12
        min_column_number = 0
        max_tries = 2
        try_number = 1
        while min_column_number != correct_number_of_output_columns and try_number <= max_tries:
            write_to_log("waiting for rapsearch server for chunk {}".format(chunk_id))
            instance_ip = wait_for_server_ip('rapsearch', key_path, remote_username, ENVIRONMENT, RAPSEARCH2_MAX_CONCURRENT, chunk_id)
            write_to_log("starting alignment for chunk %s on machine %s" % (chunk_id, instance_ip))
            execute_command_realtime_stdout(remote_command(commands, key_path, remote_username, instance_ip))
            write_to_log("finished alignment for chunk %s" % chunk_id)
            # check if every row has correct number of columns (12) in the output file on the remote machine
            verification_command = "grep -v '^#' %s" % output_path # first, remove header lines starting with '#'
            verification_command += " | awk '{print NF}' | sort -nu | head -n 1"
            min_column_number_string = execute_command_with_output(remote_command(verification_command, key_path, remote_username, instance_ip))
            min_column_number = interpret_min_column_number_string(min_column_number_string, correct_number_of_output_columns, try_number)
            try_number += 1
        # move output from remote machine to s3
        assert min_column_number == correct_number_of_output_columns, "Chunk %s output corrupt; not copying to S3. Re-start pipeline to try again." % chunk_id
        upload_command = "aws s3 cp --quiet %s %s/;" % (output_path, SAMPLE_S3_OUTPUT_CHUNKS_PATH)
        execute_command(remote_command(upload_command, key_path, remote_username, instance_ip))
        with iostream:
            execute_command(scp(key_path, remote_username, instance_ip, output_path, CHUNKS_RESULT_DIR + "/" + outfile_basename))
    return os.path.join(CHUNKS_RESULT_DIR, outfile_basename)


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


def run_rapsearch2_remotely(input_fasta, lazy_run):
    output_file = os.path.join(SAMPLE_S3_OUTPUT_PATH, RAPSEARCH2_OUT)
    key_path = fetch_key(ENVIRONMENT)
    remote_username = "ec2-user"
    remote_home_dir = "/home/%s" % remote_username
    remote_work_dir = "%s/data/batch-pipeline-workdir/%s" % (remote_home_dir, SAMPLE_NAME)
    remote_index_dir = "%s/references/nr_rapsearch" % remote_home_dir
    # split file:
    chunk_nlines = 2*RAPSEARCH_CHUNK_SIZE
    part_suffix, input_chunks = chunk_input([input_fasta], chunk_nlines, RAPSEARCH_CHUNK_SIZE)
    # process chunks:
    chunk_output_files = [None] * len(input_chunks)
    chunk_threads = []
    mutex = threading.RLock()
    iii = list(enumerate(input_chunks))
    random.shuffle(iii)
    chunks_in_flight = chunks_in_flight_rapsearch
    for n, chunk_input_files in iii:
        chunks_in_flight.acquire()
        check_for_errors(mutex, chunk_output_files, input_chunks, "rapsearch")
        t = threading.Thread(
            target=run_chunk_wrapper,
            args=[chunks_in_flight, chunk_output_files, n, mutex, run_rapsearch_chunk,
                  # Arguments passed to run_rapsearch_chunk
                  [part_suffix, remote_home_dir, remote_index_dir, remote_work_dir,
                   remote_username, chunk_input_files[0], key_path, lazy_run]]
        )
        t.start()
        chunk_threads.append(t)
    for ct in chunk_threads:
        ct.join()
        check_for_errors(mutex, chunk_output_files, input_chunks, "rapsearch")
    assert None not in chunk_output_files
    # merge output chunks:
    with iostream:
        concatenate_files(chunk_output_files, os.path.join(RESULT_DIR, RAPSEARCH2_OUT))
        execute_command("aws s3 cp --quiet %s/%s %s" % (RESULT_DIR, RAPSEARCH2_OUT, output_file))


def run_generate_taxid_outputs_from_m8(annotated_m8, taxon_counts_csv_file, taxon_counts_json_file,
                                       count_type, e_value_type, stats):

    lineage_path = fetch_reference(LINEAGE_SHELF)
    lineage_map = shelve.open(lineage_path)

    species_concordant, genus_concordant, family_concordant = generate_tax_counts_from_m8(annotated_m8, e_value_type,
                                                                                          taxon_counts_csv_file, lineage_map)
    write_to_log("generated taxon counts from m8")
    generate_json_from_taxid_counts(taxon_counts_csv_file, taxon_counts_json_file, count_type,
                                    lineage_map, species_concordant, genus_concordant, family_concordant, stats)
    write_to_log("generated JSON file from taxon counts")
    # move the output back to S3
    execute_command("aws s3 cp --quiet %s %s/" % (taxon_counts_csv_file, SAMPLE_S3_OUTPUT_PATH))
    execute_command("aws s3 cp --quiet %s %s/" % (taxon_counts_json_file, SAMPLE_S3_OUTPUT_PATH))

def run_combine_json_outputs(input_json_1, input_json_2, output_json, stats):
    combine_pipeline_output_json(input_json_1, input_json_2, output_json, stats)
    write_to_log("finished job")
    # move it the output back to S3
    execute_command("aws s3 cp --quiet %s %s/" % (output_json, SAMPLE_S3_OUTPUT_PATH))


def run_generate_unidentified_fasta(input_fa, output_fa):
    subprocess.check_output("grep -A 1 '>NR::NT::' %s | sed '/^--$/d' > %s" % (input_fa, output_fa), shell=True)
    write_to_log("finished job")
    execute_command("aws s3 cp --quiet %s %s/" % (output_fa, SAMPLE_S3_OUTPUT_PATH))


def get_alignment_version_s3_path():
    return os.path.join(SAMPLE_S3_OUTPUT_PATH, VERSION_OUT)


def run_stage2(lazy_run=True):
    # make local directories
    execute_command("mkdir -p %s %s %s %s %s" % (SAMPLE_DIR, FASTQ_DIR, RESULT_DIR, CHUNKS_RESULT_DIR, REF_DIR))

    # configure logger
    log_file = "%s/%s.%s.txt" % (RESULT_DIR, LOGS_OUT_BASENAME, AWS_BATCH_JOB_ID)
    configure_logger(log_file)

    # Download input files
    input1_s3_path = os.path.join(SAMPLE_S3_INPUT_PATH, EXTRACT_UNMAPPED_FROM_SAM_OUT1)
    input2_s3_path = os.path.join(SAMPLE_S3_INPUT_PATH, EXTRACT_UNMAPPED_FROM_SAM_OUT2)
    input3_s3_path = os.path.join(SAMPLE_S3_INPUT_PATH, EXTRACT_UNMAPPED_FROM_SAM_OUT3)
    input1_present = check_s3_file_presence(input1_s3_path)
    input2_present = check_s3_file_presence(input2_s3_path)
    input3_present = check_s3_file_presence(input3_s3_path)
    if input1_present and input2_present and input3_present:
        # output of previous stage, paired-end reads
        if SAMPLE_S3_INPUT_PATH != SAMPLE_S3_OUTPUT_PATH:
            execute_command("aws s3 cp --quiet %s %s/" % (input1_s3_path, SAMPLE_S3_OUTPUT_PATH))
            execute_command("aws s3 cp --quiet %s %s/" % (input2_s3_path, SAMPLE_S3_OUTPUT_PATH))
            execute_command("aws s3 cp --quiet %s %s/" % (input3_s3_path, SAMPLE_S3_OUTPUT_PATH))
        gsnapl_input_files = [EXTRACT_UNMAPPED_FROM_SAM_OUT1, EXTRACT_UNMAPPED_FROM_SAM_OUT2]
        before_file_name_for_log = os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_SAM_OUT1)
        before_file_type_for_log = "fasta_paired"
        # Import existing job stats
        stats = StatsFile(STATS_OUT, RESULT_DIR, SAMPLE_S3_INPUT_PATH, SAMPLE_S3_OUTPUT_PATH)
        stats.load_from_s3()
    elif input1_present and input3_present:
        # output of previous stage, non-paired-end
        if SAMPLE_S3_INPUT_PATH != SAMPLE_S3_OUTPUT_PATH:
            execute_command("aws s3 cp --quiet %s %s/" % (input1_s3_path, SAMPLE_S3_OUTPUT_PATH))
            execute_command("aws s3 cp --quiet %s %s/" % (input3_s3_path, SAMPLE_S3_OUTPUT_PATH))
        gsnapl_input_files = [EXTRACT_UNMAPPED_FROM_SAM_OUT1]
        before_file_name_for_log = os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_SAM_OUT1)
        before_file_type_for_log = "fasta_paired"
        # Import existing job stats
        stats = StatsFile(STATS_OUT, RESULT_DIR, SAMPLE_S3_INPUT_PATH, SAMPLE_S3_OUTPUT_PATH)
        stats.load_from_s3()
    else:
        # in case where previous stage was skipped, go back to original input
        command = "aws s3 ls %s/ | grep '\\.%s$'" % (SAMPLE_S3_FASTQ_PATH, FILE_TYPE)
        output = execute_command_with_output(command).rstrip().split("\n")
        for line in output:
            m = re.match(".*?([^ ]*." + re.escape(FILE_TYPE) + ")", line)
            if m:
                execute_command("aws s3 cp --quiet %s/%s %s/" % (SAMPLE_S3_FASTQ_PATH, m.group(1), FASTQ_DIR))
            else:
                write_to_log("%s doesn't match %s" % (line, FILE_TYPE))
        fastq_files = execute_command_with_output("ls %s/*.%s" % (FASTQ_DIR, FILE_TYPE)).rstrip().split("\n")
        # prepare files for gsnap
        cleaned_files, before_file_type_for_log = clean_direct_gsnapl_input(fastq_files)
        before_file_name_for_log = cleaned_files[0]
        gsnapl_input_files = [os.path.basename(f) for f in cleaned_files]
        # make combined fasta needed later
        merged_fasta = os.path.join(RESULT_DIR, EXTRACT_UNMAPPED_FROM_SAM_OUT3)
        generate_merged_fasta(cleaned_files, merged_fasta)
        execute_command("aws s3 cp --quiet %s %s/" % (merged_fasta, SAMPLE_S3_OUTPUT_PATH))
        # Create new stats
        stats = StatsFile(STATS_OUT, RESULT_DIR, None, SAMPLE_S3_OUTPUT_PATH)

    # Make sure there are no tabs in sequence names, since tabs are used as a delimiter in m8 files
    files_to_collapse_basenames = gsnapl_input_files + [EXTRACT_UNMAPPED_FROM_SAM_OUT3]
    collapsed_files = ["%s/nospace.%s" % (RESULT_DIR, f) for f in files_to_collapse_basenames]
    for file_basename in files_to_collapse_basenames:
        execute_command("aws s3 cp --quiet %s/%s %s/" % (SAMPLE_S3_OUTPUT_PATH, file_basename, RESULT_DIR))
    remove_whitespace_from_files([os.path.join(RESULT_DIR, file_basename) for file_basename in files_to_collapse_basenames],
                                 ";", collapsed_files)
    for filename in collapsed_files:
        execute_command("aws s3 cp --quiet %s %s/" % (filename, SAMPLE_S3_OUTPUT_PATH))
    gsnapl_input_files = [os.path.basename(f) for f in collapsed_files[:-1]]
    merged_fasta = collapsed_files[-1]

    # subsample if specified
    if SUBSAMPLE:
        target_n_reads = int(SUBSAMPLE)
        subsampled_gsnapl_input_files, subsampled_merged_fasta = subsample_fastas(gsnapl_input_files, os.path.basename(merged_fasta), target_n_reads)
        gsnapl_input_files = subsampled_gsnapl_input_files
        merged_fasta = os.path.join(RESULT_DIR, subsampled_merged_fasta)
        for f in gsnapl_input_files:
            execute_command("aws s3 cp --quiet %s/%s %s/" % (RESULT_DIR, f, SAMPLE_S3_OUTPUT_PATH))
        execute_command("aws s3 cp --quiet %s %s/" % (merged_fasta, SAMPLE_S3_OUTPUT_PATH))


    thread_success = {}
    thread_success_lock = threading.RLock()

    def run_gsnap_and_dependent_steps():
        # run gsnap remotely
        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "GSNAPL",
             "version_file_s3": GSNAP_VERSION_FILE_S3,
             "output_version_file": os.path.join(RESULT_DIR, VERSION_OUT)})
        run_and_log_s3(
            logparams,
            TARGET_OUTPUTS["run_gsnapl_remotely"],
            lazy_run,
            SAMPLE_S3_OUTPUT_PATH,
            run_gsnapl_remotely,
            gsnapl_input_files,
            lazy_run)
        stats.count_reads("run_gsnapl_remotely",
                          before_filename=before_file_name_for_log,
                          before_filetype=before_file_type_for_log,
                          after_filename=os.path.join(RESULT_DIR, GSNAPL_DEDUP_OUT),
                          after_filetype="m8")

        # run_annotate_gsnapl_m8_with_taxids
        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "annotate gsnapl m8 with taxids"})
        run_and_log_eager(
            logparams,
            TARGET_OUTPUTS["run_annotate_m8_with_taxids__1"],
            run_annotate_m8_with_taxids,
            os.path.join(RESULT_DIR, GSNAPL_DEDUP_OUT),
            os.path.join(RESULT_DIR, ANNOTATE_GSNAPL_M8_WITH_TAXIDS_OUT))

        # run_generate_taxid_annotated_fasta_from_m8
        logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "generate taxid annotated fasta from m8"})
        run_and_log_eager(
            logparams,
            TARGET_OUTPUTS["run_generate_taxid_annotated_fasta_from_m8__1"],
            run_generate_taxid_annotated_fasta_from_m8,
            os.path.join(RESULT_DIR, GSNAPL_DEDUP_OUT),
            merged_fasta,
            os.path.join(RESULT_DIR, GENERATE_TAXID_ANNOTATED_FASTA_FROM_M8_OUT),
            'NT')

        if SKIP_DEUTERO_FILTER:
            next_input = ANNOTATE_GSNAPL_M8_WITH_TAXIDS_OUT
        else:
            logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "filter deuterostomes from m8__1"})
            run_and_log_eager(
                logparams, TARGET_OUTPUTS["run_filter_deuterostomes_from_m8__1"],
                run_filter_deuterostomes_from_m8,
                os.path.join(RESULT_DIR, ANNOTATE_GSNAPL_M8_WITH_TAXIDS_OUT),
                os.path.join(RESULT_DIR, FILTER_DEUTEROSTOMES_FROM_NT_M8_OUT))
            stats.count_reads("run_filter_deuterostomes_from_m8__1",
                              before_filename=os.path.join(RESULT_DIR, ANNOTATE_GSNAPL_M8_WITH_TAXIDS_OUT),
                              before_filetype="m8",
                              after_filename=os.path.join(RESULT_DIR, FILTER_DEUTEROSTOMES_FROM_NT_M8_OUT),
                              after_filetype="m8")
            next_input = FILTER_DEUTEROSTOMES_FROM_NT_M8_OUT


        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "generate taxid outputs from m8"})
        run_and_log_eager(
            logparams, TARGET_OUTPUTS["run_generate_taxid_outputs_from_m8__1"],
            run_generate_taxid_outputs_from_m8,
            os.path.join(RESULT_DIR, next_input),
            os.path.join(RESULT_DIR, NT_M8_TO_TAXID_COUNTS_FILE_OUT),
            os.path.join(RESULT_DIR, NT_TAXID_COUNTS_TO_JSON_OUT),
            'NT',
            'raw',
            stats)

        with thread_success_lock:
            thread_success["gsnap"] = True

    def run_rapsearch2():
        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "RAPSearch2",
             "version_file_s3": RAPSEARCH_VERSION_FILE_S3, "output_version_file": os.path.join(RESULT_DIR, VERSION_OUT)})
        run_and_log_s3(
            logparams,
            TARGET_OUTPUTS["run_rapsearch2_remotely"],
            lazy_run,
            SAMPLE_S3_OUTPUT_PATH,
            run_rapsearch2_remotely,
            merged_fasta,
            lazy_run)
        stats.count_reads("run_rapsearch2_remotely",
                          before_filename=os.path.join(RESULT_DIR, merged_fasta),
                          before_filetype="fasta",
                          after_filename=os.path.join(RESULT_DIR, RAPSEARCH2_OUT),
                          after_filetype="m8")

        with thread_success_lock:
            thread_success["rapsearch2"] = True

    def run_additional_steps():
        # run_annotate_m8_with_taxids
        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "annotate m8 with taxids"})
        run_and_log_eager(
            logparams,
            TARGET_OUTPUTS["run_annotate_m8_with_taxids__2"],
            run_annotate_m8_with_taxids,
            os.path.join(RESULT_DIR, RAPSEARCH2_OUT),
            os.path.join(RESULT_DIR, ANNOTATE_RAPSEARCH2_M8_WITH_TAXIDS_OUT))

        if SKIP_DEUTERO_FILTER:
            next_input = ANNOTATE_RAPSEARCH2_M8_WITH_TAXIDS_OUT
        else:
            logparams = return_merged_dict(DEFAULT_LOGPARAMS, {"title": "filter deuterostomes from m8"})
            run_and_log_eager(
                logparams,
                TARGET_OUTPUTS["run_filter_deuterostomes_from_m8__2"],
                run_filter_deuterostomes_from_m8,
                os.path.join(RESULT_DIR, ANNOTATE_RAPSEARCH2_M8_WITH_TAXIDS_OUT),
                os.path.join(RESULT_DIR, FILTER_DEUTEROSTOMES_FROM_NR_M8_OUT))
            stats.count_reads("run_filter_deuterostomes_from_m8__2",
                              before_filename=os.path.join(RESULT_DIR, ANNOTATE_RAPSEARCH2_M8_WITH_TAXIDS_OUT),
                              before_filetype="fasta",
                              after_filename=os.path.join(RESULT_DIR, FILTER_DEUTEROSTOMES_FROM_NR_M8_OUT),
                              after_filetype="m8")
            next_input = FILTER_DEUTEROSTOMES_FROM_NR_M8_OUT

        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "generate taxid outputs from m8"})
        run_and_log_eager(
            logparams,
            TARGET_OUTPUTS["run_generate_taxid_outputs_from_m8__2"],
            run_generate_taxid_outputs_from_m8,
            os.path.join(RESULT_DIR, next_input),
            os.path.join(RESULT_DIR, NR_M8_TO_TAXID_COUNTS_FILE_OUT),
            os.path.join(RESULT_DIR, NR_TAXID_COUNTS_TO_JSON_OUT),
            'NR',
            'log10',
            stats)

        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "combine JSON outputs"})
        run_and_log_eager(
            logparams,
            TARGET_OUTPUTS["run_combine_json_outputs"],
            run_combine_json_outputs,
            RESULT_DIR + '/' + NT_TAXID_COUNTS_TO_JSON_OUT,
            RESULT_DIR + '/' + NR_TAXID_COUNTS_TO_JSON_OUT,
            RESULT_DIR + '/' + COMBINED_JSON_OUT,
            stats)

        with thread_success_lock:
            thread_success["additional_steps"] = True

    def run_annotation_steps():
        # run_generate_taxid_annotated_fasta_from_m8
        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "generate taxid annotated fasta from m8"})
        run_and_log_eager(
            logparams,
            TARGET_OUTPUTS["run_generate_taxid_annotated_fasta_from_m8__2"],
            run_generate_taxid_annotated_fasta_from_m8,
            RESULT_DIR + '/' + RAPSEARCH2_OUT,
            RESULT_DIR + '/' + GENERATE_TAXID_ANNOTATED_FASTA_FROM_M8_OUT,
            RESULT_DIR + '/' + GENERATE_TAXID_ANNOTATED_FASTA_FROM_RAPSEARCH2_M8_OUT,
            'NR')

        logparams = return_merged_dict(
            DEFAULT_LOGPARAMS,
            {"title": "generate FASTA of unidentified reads"})
        run_and_log_eager(
            logparams,
            TARGET_OUTPUTS["run_generate_unidentified_fasta"],
            run_generate_unidentified_fasta,
            RESULT_DIR + '/' + GENERATE_TAXID_ANNOTATED_FASTA_FROM_RAPSEARCH2_M8_OUT,
            RESULT_DIR + '/' + UNIDENTIFIED_FASTA_OUT)
        stats.count_reads("run_generate_unidentified_fasta",
                          before_filename=os.path.join(RESULT_DIR, GENERATE_TAXID_ANNOTATED_FASTA_FROM_RAPSEARCH2_M8_OUT),
                          before_filetype="fasta",
                          after_filename=os.path.join(RESULT_DIR, UNIDENTIFIED_FASTA_OUT),
                          after_filetype="fasta")


        with thread_success_lock:
            thread_success["annotation"] = True

    t_gsnap = threading.Thread(target=run_gsnap_and_dependent_steps)
    t_rapsearch2 = threading.Thread(target=run_rapsearch2)
    t_additional_steps = threading.Thread(target=run_additional_steps)
    t_annotation = threading.Thread(target=run_annotation_steps)

    def thread_succeded(name):
        with thread_success_lock:
            return thread_success.get(name)

    t_gsnap.start()
    t_rapsearch2.start()
    t_gsnap.join()
    assert thread_succeded("gsnap")
    t_rapsearch2.join()
    assert thread_succeded("rapsearch2")

    t_additional_steps.start()
    t_annotation.start()
    t_additional_steps.join()
    assert thread_succeded("additional_steps")
    t_annotation.join()
    assert thread_succeded("annotation")

    # copy log file -- after work is done
    stats.save_to_s3()
    upload_log_file(SAMPLE_S3_OUTPUT_PATH)
