import os
import subprocess
import json
import shelve
import logging
import idseq_pipeline.commands.accessionid2seq_functions as accessionid2seq_functions
from .common import *  #pylint: disable=wildcard-import

# data directories
# from common import ROOT_DIR
# from common import REF_DIR
DEST_DIR = ROOT_DIR + '/idseq/data'  # generated data go here
TEMP_DIR = ROOT_DIR + '/tmp'  # tmp directory with a lot of space for sorting large files

# arguments from environment variables
INPUT_BUCKET = os.environ.get('INPUT_BUCKET')
OUTPUT_BUCKET = os.environ.get('OUTPUT_BUCKET')
AWS_BATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID', 'local')
SAMPLE_S3_INPUT_PATH = INPUT_BUCKET.rstrip('/')
SAMPLE_S3_OUTPUT_PATH = OUTPUT_BUCKET.rstrip('/')
sample_name = SAMPLE_S3_INPUT_PATH[5:].rstrip('/').replace('/', '-')
SAMPLE_DIR = DEST_DIR + '/' + sample_name
INPUT_DIR = SAMPLE_DIR + '/inputs'
RESULT_DIR = SAMPLE_DIR + '/results'
DEFAULT_LOG_PARAMS = {"sample_s3_output_path": SAMPLE_S3_OUTPUT_PATH}

NT_LOC_DB = os.environ.get('NT_LOC_DB',
                           "s3://idseq-database/20170824/blast_db/nt_loc.db")
NT_DB = os.environ.get('NT_DB', "s3://idseq-database/20170824/blast_db/nt")

# input files
ACCESSION_ANNOTATED_FASTA = 'accessions.rapsearch2.gsnapl.fasta'
GSNAP_M8_FILE = 'dedup.multihit.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8'
SUMMARY_MULTIHIT_GSNAPL_OUT = 'summary.multihit.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.tab'
SUMMARY_MULTIHIT_RAPSEARCH_OUT = 'summary.multihit.rapsearch2.filter.deuterostomes.taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.tab'

# output files
TAXID_ANNOT_FASTA = 'taxid_annot.fasta'
TAXID_ANNOT_SORTED_FASTA_NT = 'taxid_annot_sorted_nt.fasta'
TAXID_ANNOT_SORTED_FASTA_NR = 'taxid_annot_sorted_nr.fasta'
TAXID_ANNOT_SORTED_FASTA_GENUS_NT = 'taxid_annot_sorted_genus_nt.fasta'
TAXID_ANNOT_SORTED_FASTA_GENUS_NR = 'taxid_annot_sorted_genus_nr.fasta'
TAXID_ANNOT_SORTED_FASTA_FAMILY_NT = 'taxid_annot_sorted_family_nt.fasta'
TAXID_ANNOT_SORTED_FASTA_FAMILY_NR = 'taxid_annot_sorted_family_nr.fasta'
TAXID_LOCATIONS_JSON_NT = 'taxid_locations_nt.json'
TAXID_LOCATIONS_JSON_NR = 'taxid_locations_nr.json'
TAXID_LOCATIONS_JSON_GENUS_NT = 'taxid_locations_genus_nt.json'
TAXID_LOCATIONS_JSON_GENUS_NR = 'taxid_locations_genus_nr.json'
TAXID_LOCATIONS_JSON_FAMILY_NT = 'taxid_locations_family_nt.json'
TAXID_LOCATIONS_JSON_FAMILY_NR = 'taxid_locations_family_nr.json'
TAXID_LOCATIONS_JSON_ALL = 'taxid_locations_combined.json'
LOGS_OUT_BASENAME = 'postprocess-log'
ALIGN_VIZ_DIR = 'align_viz'

# target outputs by task
TARGET_OUTPUTS = {
    "run_generate_taxid_fasta_from_hit_summaries":
    [os.path.join(RESULT_DIR, TAXID_ANNOT_FASTA)],
    "run_generate_taxid_locator__1": [
        os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_NT),
        os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_NT)
    ],
    "run_generate_taxid_locator__2": [
        os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_NR),
        os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_NR)
    ],
    "run_generate_taxid_locator__3": [
        os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_GENUS_NT),
        os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_GENUS_NT)
    ],
    "run_generate_taxid_locator__4": [
        os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_GENUS_NR),
        os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_GENUS_NR)
    ],
    "run_generate_taxid_locator__5": [
        os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_FAMILY_NT),
        os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_FAMILY_NT)
    ],
    "run_generate_taxid_locator__6": [
        os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_FAMILY_NR),
        os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_FAMILY_NR)
    ],
    "run_combine_json": [os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_ALL)],
    "run_generate_align_viz":
    [os.path.join(RESULT_DIR, "%s.summary" % ALIGN_VIZ_DIR)]
}

# references
# from common import ACCESSION2TAXID


# processing functions
def remove_annotation(read_id):
    result = re.sub(r'NT:(.*?):', '', read_id)
    result = re.sub(r'NR:(.*?):', '', result)
    return result


def parse_hits(hit_summary_files):
    # Return map of {NT, NR} => read_id => (hit_taxid_str, hit_level_str)
    valid_hits = {}
    for count_type, summary_file in hit_summary_files.items():
        hits = {}
        with open(summary_file, "rb") as sf:
            for hit_line in sf:
                hit_line_columns = hit_line.strip().split("\t")
                if len(hit_line_columns) >= 3:
                    hit_read_id = hit_line_columns[0]
                    hit_level_str = hit_line_columns[1]
                    hit_taxid_str = hit_line_columns[2]
                    hits[hit_read_id] = (hit_taxid_str, hit_level_str)
        valid_hits[count_type] = hits
    return valid_hits


def generate_taxid_fasta_from_hit_summaries(
        input_fasta_file, hit_summary_files, lineagePath, output_fasta_file):
    lineage_map = shelve.open(lineagePath)
    valid_hits = parse_hits(hit_summary_files)

    def get_valid_lineage(read_id, count_type):
        # If the read aligned to something, then it would be present in the summary file
        # for count type, and correspondingly in valid_hits[count_type], even if the
        # hits disagree so much that the "valid_hits" entry is just ("-1", "-1").
        # If the read didn't align to anything, we also represent that with ("-1", "-1).
        # This ("-1", "-1) gets translated to NULL_LINEAGE.
        hit_taxid_str, hit_level_str = valid_hits[count_type].get(
            read_id, ("-1", "-1"))
        hit_lineage = lineage_map.get(hit_taxid_str, NULL_LINEAGE)
        return validate_taxid_lineage(hit_lineage, hit_taxid_str,
                                      hit_level_str)

    input_fasta_f = open(input_fasta_file, 'rb')
    output_fasta_f = open(output_fasta_file, 'wb')
    sequence_name = input_fasta_f.readline()
    sequence_data = input_fasta_f.readline()
    while len(sequence_name) > 0 and len(sequence_data) > 0:
        accession_annotated_read_id = sequence_name.rstrip().lstrip(
            '>'
        )  # example read_id: "NR::NT:CP010376.2:NB501961:14:HM7TLBGX2:1:23109:12720:8743/2"
        read_id = accession_annotated_read_id.split(":", 4)[-1]

        nr_taxid_species, nr_taxid_genus, nr_taxid_family = get_valid_lineage(
            read_id, 'NR')
        nt_taxid_species, nt_taxid_genus, nt_taxid_family = get_valid_lineage(
            read_id, 'NT')

        new_read_name = (
            'family_nr:' + nr_taxid_family + ':family_nt:' + nt_taxid_family +
            ':genus_nr:' + nr_taxid_genus + ':genus_nt:' + nt_taxid_genus +
            ':species_nr:' + nr_taxid_species + ':species_nt:' +
            nt_taxid_species + ':' + accession_annotated_read_id)
        output_fasta_f.write(">%s\n" % new_read_name)
        output_fasta_f.write(sequence_data)
        sequence_name = input_fasta_f.readline()
        sequence_data = input_fasta_f.readline()
    input_fasta_f.close()
    output_fasta_f.close()


def get_taxid(sequence_name, taxid_field):
    parts = sequence_name.replace('>', ':').split(":%s:" % taxid_field)
    if len(parts) <= 1:
        # sequence_name empty or taxid_field not found
        return 'none'
    taxid = parts[1].split(":")[0]
    # example sequence_name: ">nr:-100:nt:684552:NR::NT:LT629734.1:HWI-ST640:828:H917FADXX:2:1101:1424:15119/1"
    return taxid


def get_taxid_field_num(taxid_field, input_fasta):
    with open(input_fasta) as f:
        sequence_name = f.readline()
    return sequence_name.replace('>', ':').split(":").index(taxid_field) + 1


def generate_taxid_locator(input_fasta, taxid_field, hit_type, output_fasta,
                           output_json):
    taxid_field_num = get_taxid_field_num(taxid_field, input_fasta)
    # put every 2-line fasta record on a single line with delimiter ":lineseparator:":
    command = "awk 'NR % 2 == 1 { o=$0 ; next } { print o \":lineseparator:\" $0 }' " + input_fasta
    # sort the records based on the field containing the taxids:
    command += " | sort -T %s --key %s --field-separator ':' --numeric-sort" % (
        TEMP_DIR, taxid_field_num)
    # split every record back over 2 lines:
    command += " | sed 's/:lineseparator:/\\n/g' > %s" % output_fasta
    subprocess.check_output(command, shell=True)
    # make json giving byte range of file corresponding to each taxid:
    taxon_sequence_locations = []
    f = open(output_fasta, 'rb')
    sequence_name = f.readline()
    sequence_data = f.readline()
    taxid = get_taxid(sequence_name, taxid_field)
    first_byte = 0
    end_byte = first_byte + len(sequence_name) + len(sequence_data)
    while len(sequence_name) > 0 and len(sequence_data) > 0:
        sequence_name = f.readline()
        sequence_data = f.readline()
        new_taxid = get_taxid(sequence_name, taxid_field)
        if new_taxid != taxid:
            # Note on boundary condition: when end of file is reached, then
            # sequence_name == '' => new_taxid=='none' => new_taxid != taxid
            # so last record will be written to output correctly.
            taxon_sequence_locations.append({
                'taxid': int(taxid),
                'first_byte': first_byte,
                'last_byte': end_byte - 1,
                'hit_type': hit_type
            })
            taxid = new_taxid
            first_byte = end_byte
            end_byte = first_byte + len(sequence_name) + len(sequence_data)
        else:
            end_byte += len(sequence_name) + len(sequence_data)
    f.close()
    with open(output_json, 'wb') as f:
        json.dump(taxon_sequence_locations, f)


def combine_json(input_json_list, output_json):
    output = []
    for input_json in input_json_list:
        with open(input_json) as f:
            output.extend(json.load(f))
    with open(output_json, 'wb') as outf:
        json.dump(output, outf)


# job functions


def run_generate_align_viz(input_fasta, input_m8, output_dir):
    nt_loc_db = fetch_reference(NT_LOC_DB)
    summary_file_name = accessionid2seq_functions.generate_alignment_viz_json(
        NT_DB, nt_loc_db, "NT", input_m8, input_fasta, output_dir)
    # copy the data over
    execute_command("aws s3 cp --quiet %s %s/align_viz --recursive" %
                    (output_dir, SAMPLE_S3_OUTPUT_PATH))
    execute_command("aws s3 cp --quiet %s %s/" % (summary_file_name,
                                                  SAMPLE_S3_OUTPUT_PATH))


def run_generate_taxid_fasta_from_hit_summaries(input_fasta, hit_summary_files,
                                                output_fasta):
    lineage_path = fetch_reference(LINEAGE_SHELF)
    generate_taxid_fasta_from_hit_summaries(input_fasta, hit_summary_files,
                                            lineage_path, output_fasta)
    logging.getLogger().info("finished job")
    execute_command("aws s3 cp --quiet %s %s/" % (output_fasta,
                                                  SAMPLE_S3_OUTPUT_PATH))


def run_generate_taxid_locator(input_fasta, taxid_field, hit_type,
                               output_fasta, output_json):
    generate_taxid_locator(input_fasta, taxid_field, hit_type, output_fasta,
                           output_json)
    logging.getLogger().info("finished job")
    execute_command("aws s3 cp --quiet %s %s/" % (output_fasta,
                                                  SAMPLE_S3_OUTPUT_PATH))
    execute_command("aws s3 cp --quiet %s %s/" % (output_json,
                                                  SAMPLE_S3_OUTPUT_PATH))


def run_combine_json(input_json_list, output_json):
    combine_json(input_json_list, output_json)
    logging.getLogger().info("finished job")
    execute_command("aws s3 cp --quiet %s %s/" % (output_json,
                                                  SAMPLE_S3_OUTPUT_PATH))


def run_stage3(lazy_run=False):
    assert not lazy_run, "run_stage3 was called with lazy_run"

    # Make data directories
    execute_command("mkdir -p %s %s %s %s" % (SAMPLE_DIR, RESULT_DIR, REF_DIR,
                                              TEMP_DIR))

    # Configure logger
    log_file = "%s/%s.%s.txt" % (RESULT_DIR, LOGS_OUT_BASENAME,
                                 AWS_BATCH_JOB_ID)
    configure_logger(log_file)

    print("Starting stage...")

    # Download input
    execute_command("aws s3 cp --quiet %s/%s %s/" %
                    (SAMPLE_S3_INPUT_PATH, ACCESSION_ANNOTATED_FASTA,
                     INPUT_DIR))
    input_file = os.path.join(INPUT_DIR, ACCESSION_ANNOTATED_FASTA)

    # Download m8
    execute_command("aws s3 cp --quiet %s/%s %s/" % (SAMPLE_S3_INPUT_PATH,
                                                     GSNAP_M8_FILE, INPUT_DIR))
    input_m8 = os.path.join(INPUT_DIR, GSNAP_M8_FILE)

    # Download hit level files
    hit_summary_files = {
        'NT': os.path.join(INPUT_DIR, SUMMARY_MULTIHIT_GSNAPL_OUT),
        'NR': os.path.join(INPUT_DIR, SUMMARY_MULTIHIT_RAPSEARCH_OUT)
    }
    for local_file in hit_summary_files.values():
        execute_command("aws s3 cp --quiet %s/%s %s/" %
                        (SAMPLE_S3_INPUT_PATH, os.path.basename(local_file),
                         INPUT_DIR))

    def default_log_param_plus_title(title):
        return return_merged_dict(DEFAULT_LOG_PARAMS, {"title": title})

    # Ex: run_and_log(log_params, target_outputs, lazy_run, func_name, *args)
    # TODO: Get rid of run_and_log pattern for simplification

    # Generate taxid fasta
    log_params = default_log_param_plus_title(
        "run_generate_taxid_fasta_from_hit_summaries")
    run_and_log(log_params,
                TARGET_OUTPUTS["run_generate_taxid_fasta_from_hit_summaries"],
                False, run_generate_taxid_fasta_from_hit_summaries, input_file,
                hit_summary_files, os.path.join(RESULT_DIR, TAXID_ANNOT_FASTA))

    # SPECIES level

    # Generate taxid locator for NT
    log_params = default_log_param_plus_title("run_generate_taxid_locator for NT")
    run_and_log(log_params, TARGET_OUTPUTS["run_generate_taxid_locator__1"],
                False, run_generate_taxid_locator,
                os.path.join(RESULT_DIR,
                             TAXID_ANNOT_FASTA), 'species_nt', 'NT',
                os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_NT),
                os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_NT))

    # Generate taxid locator for NR
    log_params = default_log_param_plus_title("run_generate_taxid_locator for NR")
    run_and_log(log_params, TARGET_OUTPUTS["run_generate_taxid_locator__2"],
                False, run_generate_taxid_locator,
                os.path.join(RESULT_DIR,
                             TAXID_ANNOT_FASTA), 'species_nr', 'NR',
                os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_NR),
                os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_NR))

    # GENUS level

    # Generate taxid locator for NT
    log_params = default_log_param_plus_title("run_generate_taxid_locator for NT")
    run_and_log(log_params, TARGET_OUTPUTS["run_generate_taxid_locator__3"],
                False, run_generate_taxid_locator,
                os.path.join(RESULT_DIR, TAXID_ANNOT_FASTA), 'genus_nt', 'NT',
                os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_GENUS_NT),
                os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_GENUS_NT))

    # Generate taxid locator for NR
    log_params = default_log_param_plus_title("run_generate_taxid_locator for NR")
    run_and_log(log_params, TARGET_OUTPUTS["run_generate_taxid_locator__4"],
                False, run_generate_taxid_locator,
                os.path.join(RESULT_DIR, TAXID_ANNOT_FASTA), 'genus_nr', 'NR',
                os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_GENUS_NR),
                os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_GENUS_NR))

    # FAMILY level

    # Generate taxid locator for NT
    log_params = default_log_param_plus_title("run_generate_taxid_locator for NT")
    run_and_log(log_params, TARGET_OUTPUTS["run_generate_taxid_locator__5"],
                False, run_generate_taxid_locator,
                os.path.join(RESULT_DIR, TAXID_ANNOT_FASTA), 'family_nt', 'NT',
                os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_FAMILY_NT),
                os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_FAMILY_NT))

    # Generate taxid locator for NR
    log_params = default_log_param_plus_title("run_generate_taxid_locator for NR")
    run_and_log(log_params, TARGET_OUTPUTS["run_generate_taxid_locator__6"],
                False, run_generate_taxid_locator,
                os.path.join(RESULT_DIR, TAXID_ANNOT_FASTA), 'family_nr', 'NR',
                os.path.join(RESULT_DIR, TAXID_ANNOT_SORTED_FASTA_FAMILY_NR),
                os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_FAMILY_NR))

    # Generate alignment visualization
    log_params = default_log_param_plus_title("run_generate_align_viz")
    run_and_log(log_params, TARGET_OUTPUTS["run_generate_align_viz"], False,
                run_generate_align_viz,
                os.path.join(RESULT_DIR,
                             TAXID_ANNOT_SORTED_FASTA_NT), input_m8,
                os.path.join(RESULT_DIR, ALIGN_VIZ_DIR))

    # Combine results
    log_params = default_log_param_plus_title("run_combine_json")
    input_files_basenames = [
        TAXID_LOCATIONS_JSON_NT, TAXID_LOCATIONS_JSON_NR,
        TAXID_LOCATIONS_JSON_GENUS_NT, TAXID_LOCATIONS_JSON_GENUS_NR,
        TAXID_LOCATIONS_JSON_FAMILY_NT, TAXID_LOCATIONS_JSON_FAMILY_NR
    ]
    input_files = [os.path.join(RESULT_DIR, f) for f in input_files_basenames]
    run_and_log(log_params, TARGET_OUTPUTS["run_combine_json"], False,
                run_combine_json, input_files,
                os.path.join(RESULT_DIR, TAXID_LOCATIONS_JSON_ALL))
