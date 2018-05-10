import os
import subprocess
import json
from .common import * #pylint: disable=wildcard-import

# data directories
# from common import ROOT_DIR
# from common import REF_DIR
DEST_DIR = ROOT_DIR + '/idseq/data' # generated data go here
TEMP_DIR = ROOT_DIR + '/tmp' # tmp directory with a lot of space for sorting large files
SPADES_DIR = ROOT_DIR + '/spades' # outputs of SPAdes assemblies go here

# parameters
ASSEMBLY_READ_THRESHOLD = 100
MAX_TAXIDS_TO_ASSEMBLE = 100
MAX_NONHOST_READS = 10**5

# arguments from environment variables
ALIGNMENT_S3_PATH = os.environ.get('ALIGNMENT_S3_PATH').rstrip('/')
POSTPROCESS_S3_PATH = os.environ.get('POSTPROCESS_S3_PATH').rstrip('/')
AWS_BATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID', 'local')
sample_name = ALIGNMENT_S3_PATH[5:].rstrip('/').replace('/', '-')
SAMPLE_DIR = DEST_DIR + '/' + sample_name
INPUT_DIR = SAMPLE_DIR + '/inputs'
RESULT_DIR = SAMPLE_DIR + '/results'

# input files
TAXON_COUNTS = 'multihit_idseq_web_sample.json'
TAXID_ANNOT_FASTA = 'taxid_annot.fasta'

# outputs
SAMPLE_S3_OUTPUT_PATH = POSTPROCESS_S3_PATH
ASSEMBLY_DIR = 'assembly'
STATUS_FILE = 'job-complete'
ASSEMBLY_LOGFILE = 'assembly.log'

def run_stage4():

    # make data directories
    execute_command("mkdir -p %s %s %s %s %s" % (SAMPLE_DIR, RESULT_DIR, REF_DIR, TEMP_DIR, SPADES_DIR))

    # download input
    execute_command("aws s3 cp --quiet %s/%s %s/" % (ALIGNMENT_S3_PATH, TAXON_COUNTS, INPUT_DIR))
    pipeline_output_json = os.path.join(INPUT_DIR, TAXON_COUNTS)
    execute_command("aws s3 cp --quiet %s/%s %s/" % (POSTPROCESS_S3_PATH, TAXID_ANNOT_FASTA, INPUT_DIR))
    full_fasta = os.path.join(INPUT_DIR, TAXID_ANNOT_FASTA)

    # run assembly
    def get_taxids_to_assemble(taxon_counts, top=MAX_TAXIDS_TO_ASSEMBLE, read_threshold=ASSEMBLY_READ_THRESHOLD):
        total_counts = {} # sum of NT and NR counts for each taxid
        filtered_total_counts = {} # taxids that satisfy count threshold
        for item in taxon_counts:
            taxid = item['tax_id']
            total_counts[taxid] = total_counts.get(taxid, 0) + item['count']
            if read_threshold is None or total_counts[taxid] >= read_threshold:
                filtered_total_counts[taxid] = total_counts[taxid]
        sorted_filtered_taxids = sorted(filtered_total_counts, key=filtered_total_counts.get, reverse=True)        
        if top is not None:
            sorted_filtered_taxids = sorted_filtered_taxids[:top]
        return sorted_filtered_taxids

    def s3_result_exists(s3_path):
        try:
            return int(execute_command_with_output("aws s3 ls %s | wc -l" % s3_path).rstrip())
        except:
            return 0

    def make_inputs_for_assembly(lazy_run=True):
        with open(pipeline_output_json) as f:
            pipeline_output = json.load(f)
        taxon_counts = pipeline_output['pipeline_output']['taxon_counts_attributes']
        sorted_taxids_to_assemble = get_taxids_to_assemble(taxon_counts)
        # Get reads for the taxids chosen for assembly
        output = {}
        hit_delimiters = ['_nt', '_nr'] # annotations are like e.g. "genus_nt:543:"
        for taxid in sorted_taxids_to_assemble:
            if lazy_run and s3_result_exists(os.path.join(SAMPLE_S3_OUTPUT_PATH, ASSEMBLY_DIR, taxid)):
                continue
            pattern = '\|'.join(['%s:%s:' % (delimiter, taxid) for delimiter in hit_delimiters])
            partial_fasta =  os.path.join(RESULT_DIR, taxid + ".fasta")
            try:
                subprocess.check_call("grep -A 1 --no-group-separator '%s' %s > %s" % (pattern, full_fasta, partial_fasta), shell=True)
                output[taxid] = partial_fasta
            except:
                print "WARNING: taxid %s was not found in the annotated fasta" % taxid
        # Also include the full fasta as an input to assembly if it is not too large
        nonhost_reads = pipeline_output['pipeline_output']['remaining_reads']
        if nonhost_reads <= MAX_NONHOST_READS:
            if not lazy_run or not s3_result_exists(os.path.join(SAMPLE_S3_OUTPUT_PATH, ASSEMBLY_DIR, 'all')):
                output['all'] = full_fasta
                sorted_taxids_to_assemble.append('all')
        sorted_taxids_to_assemble = [item for item in sorted_taxids_to_assemble if item in output.keys()]
        return output, sorted_taxids_to_assemble

    def length_without_newlines(sequence):
        return len(sequence.replace("\n",""))

    def clean_scaffolds(spades_scaffolds_fasta, min_contig_length, cleaned_fasta):
        def conditional_write_record(name, sequence, output_file):
            if length_without_newlines(sequence) > min_contig_length:
                output_file.write(name)
                output_file.write(sequence)
        with open(spades_scaffolds_fasta, 'rb') as input_f:
            with open(cleaned_fasta, 'wb') as output_f:
                name = ''
                sequence = ''
                line = input_f.readline()
                while len(line) > 0:
                   if line.startswith('>'):
                       conditional_write_record(name, sequence, output_f)
                       name = line
                       sequence = ''
                   else:
                      sequence += line
                   line = input_f.readline()
                conditional_write_record(name, sequence, output_f)

    def max_read_length(input_fasta, n=2):
        ''' get read length by looking at first n records '''
        record_number = 0
        read_length = 0
        max_length = 0
        with open(input_fasta, 'rb') as f:
            while record_number <= n:
                line = f.readline()
                if line.startswith('>'):
                    max_length = max(max_length, read_length)
                    read_length = 0
                    i += 1
                else:
                   read_length += length_without_newlines(line)
        return max_length

    def spades(input_fasta, output_fasta):
        tmp_output_dir = input_fasta + "_temp_output"
        tmp_output_file = tmp_output_dir + "/scaffolds.fasta"
        try:
            execute_command("spades.py -s %s -o %s -m 60 -t 32 --only-assembler" % (input_fasta, tmp_output_dir))
            subprocess.check_call("mv %s %s" % (tmp_output_file, output_fasta), shell=True)
            return True
        except:
            traceback.print_exc()
            return False

    inputs, sorted_taxids = make_inputs_for_assembly()
    assembly_logfile = os.path.join(SAMPLE_S3_OUTPUT_PATH, ASSEMBLY_LOGFILE)
    for taxid in sorted_taxids:
        input_fasta = inputs[taxid]
        spades_output = os.path.join(RESULT_DIR, taxid + ".scaffolds.fasta")
        output_fasta = os.path.join(RESULT_DIR, taxid + ".cleaned-scaffolds.fasta")
        if spades(input_fasta, spades_output):
            clean_scaffolds(spades_output, max_read_length(input_fasta), output_fasta)
            output_s3 = "%s/%s/%s" % (SAMPLE_S3_OUTPUT_PATH, ASSEMBLY_DIR, taxid)
            execute_command("aws s3 cp --quiet %s %s" % (output_fasta, output_s3))

    # Finally, upload status file so web app knows we're done
    execute_command("echo '' | aws s3 cp --quiet - %s/%s-%s" % (SAMPLE_S3_OUTPUT_PATH, ASSEMBLY_DIR, STATUS_FILE))
