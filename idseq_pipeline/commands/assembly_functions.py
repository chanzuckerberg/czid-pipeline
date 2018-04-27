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

# arguments from environment variables
ALIGNMENT_S3_PATH = os.environ.get('ALIGNMENT_S3_PATH').rstrip('/')
POSTPROCESS_S3_PATH = os.environ.get('POSTPROCESS_S3_PATH').rstrip('/')
AWS_BATCH_JOB_ID = os.environ.get('AWS_BATCH_JOB_ID', 'local')
sample_name = ALIGNMENT_S3_PATH[5:].rstrip('/').replace('/', '-')
SAMPLE_DIR = DEST_DIR + '/' + sample_name
INPUT_DIR = SAMPLE_DIR + '/inputs'
RESULT_DIR = SAMPLE_DIR + '/results'

# input files
TAXON_COUNTS = 'idseq_web_sample.json'
TAXID_ANNOT_FASTA = 'taxid_annot.fasta'

# outputs
SAMPLE_S3_OUTPUT_PATH = POSTPROCESS_S3_PATH
ASSEMBLY_DIR = 'assembly'
STATUS_FILE = 'job-complete'

def run_stage4():

    # make data directories
    execute_command("mkdir -p %s %s %s %s %s" % (SAMPLE_DIR, RESULT_DIR, REF_DIR, TEMP_DIR, SPADES_DIR))

    # download input
    execute_command("aws s3 cp --quiet %s/%s %s/" % (ALIGNMENT_S3_PATH, TAXON_COUNTS, INPUT_DIR))
    pipeline_output_json = os.path.join(INPUT_DIR, TAXON_COUNTS)
    execute_command("aws s3 cp --quiet %s/%s %s/" % (POSTPROCESS_S3_PATH, TAXID_ANNOT_FASTA, INPUT_DIR))
    full_fasta = os.path.join(INPUT_DIR, TAXID_ANNOT_FASTA)

    # run assembly
    def get_taxids_to_assemble(taxon_counts):
        total_counts = {} # sum of NT and NR counts for each taxid
        filtered_total_counts = {} # taxids that satisfy count threshold
        for item in taxon_counts:
            taxid = item['tax_id']
            total_counts[taxid] = total_counts.get(taxid, 0) + item['count']
            if total_counts[taxid] >= ASSEMBLY_READ_THRESHOLD:
                filtered_total_counts[taxid] = total_counts[taxid]
        sorted_filtered_taxids = sorted(filtered_total_counts, key=filtered_total_counts.get, reverse=True)        
        return sorted_filtered_taxids

    def make_inputs_for_assembly():
        # Get taxids to assemble based on criteria from report, currently: taxids with read count >= ASSEMBLY_READ_THRESHOLD
        with open(pipeline_output_json) as f:
            pipeline_output = json.load(f)
        taxon_counts = pipeline_output['pipeline_output']['taxon_counts_attributes']
        sorted_taxids_to_assemble = get_taxids_to_assemble(taxon_counts)
        # Get reads for those taxids
        output = {}
        hit_delimiters = ['_nt', '_nr'] # annotations are like e.g. "genus_nt:543:"
        for taxid in sorted_taxids_to_assemble:
            pattern = '\|'.join(['%s:%s:' % (delimiter, taxid) for delimiter in hit_delimiters])
            partial_fasta =  os.path.join(RESULT_DIR, taxid + ".fasta")
            try:
                subprocess.check_call("grep -A 1 --no-group-separator '%s' %s > %s" % (pattern, full_fasta, partial_fasta), shell=True)
                output[taxid] = partial_fasta
            except:
                print "WARNING: taxid %s was not found in the annotated fasta" % taxid
        # Also include the full fasta as an input to assembly
        output['all'] = full_fasta
        sorted_taxids_to_assemble.append('all')
        sorted_taxids_to_assemble = [item in sorted_taxids_to_assemble if item in output.keys()]
        return output, sorted_taxids_to_assemble

    def length_without_newlines(sequence):
        return len(sequence.replace("\n",""))

    def read_fasta(input_fasta):
        ''' Read fasta file into dictionary. Newlines remain part of the keys/values. '''
        read_name = None
        d = {}
        read_order = []
        with open(input_fasta, 'rb') as f:
            for line in f:
                if line.startswith('>'):
                    read_name = line
                    read_order += [read_name]
                else:
                    d[read_name] = d.get(read_name, '') + line
        return d, read_order

    def write_fasta(fasta_dict, read_order, output_file):
        with open(output_file, 'wb') as f:
            for read_name in read_order:
                f.write(read_name + fasta_dict[read_name])

    def clean_scaffolds(spades_scaffolds_fasta, min_contig_length, cleaned_fasta):
        fasta_dict, read_order = read_fasta(spades_scaffolds_fasta)
        for name, sequence in fasta_dict.iteritems():
            if length_without_newlines(sequence) <= min_contig_length:
                read_order.remove(name)
        write_fasta(fasta_dict, read_order, cleaned_fasta)

    def max_read_length(input_fasta):
        fasta_dict, _read_order = read_fasta(input_fasta)
        return max([length_without_newlines(sequence) for sequence in fasta_dict.values()])

    def spades(input_fasta, output_fasta):
        tmp_output_dir = input_fasta + "_temp_output"
        tmp_output_file = tmp_output_dir + "/scaffolds.fasta"
        try:
            execute_command_realtime_stdout("spades.py -s %s -o %s -m 60 -t 32 --only-assembler" % (input_fasta, tmp_output_dir))
            subprocess.check_call("mv %s %s" % (tmp_output_file, output_fasta), shell=True)
            return True
        except:
            return False

    inputs, sorted_taxids = make_inputs_for_assembly()
    for taxid in sorted_taxids:
        input_fasta = inputs[taxid]
        spades_output = os.path.join(RESULT_DIR, taxid + ".scaffolds.fasta")
        output_fasta = os.path.join(RESULT_DIR, taxid + ".cleaned-scaffolds.fasta")
        if spades(input_fasta, spades_output):
            clean_scaffolds(spades_output, max_read_length(input_fasta), output_fasta)
            if taxid == 'all':
                output_s3 = "%s/%s-%s" % (SAMPLE_S3_OUTPUT_PATH, ASSEMBLY_DIR, taxid)
            else:
                output_s3 = "%s/%s/%s" % (SAMPLE_S3_OUTPUT_PATH, ASSEMBLY_DIR, taxid)
            execute_command("aws s3 cp --quiet %s %s" % (output_fasta, output_s3))

    # Finally, upload status file so web app knows we're done
    execute_command("echo '' | aws s3 cp --quiet - %s/%s-%s" % (SAMPLE_S3_OUTPUT_PATH, ASSEMBLY_DIR, STATUS_FILE))
