"""
Script file for generating reference updates
"""

from idseq_pipeline.commands.common import *

local_work_dir = "idseq_pipeline_temp"
ncbitool_s3_path = "s3://czbiohub-infectious-disease/ncbitool"


def main():
    ##### PARAMETERS TO EDIT #####

    ## To pull NCBI references through ncbitool, which archives and records versions, set URL_PREFIX=''.
    # To pull the latest files from NCBI without using ncbitool (no version information), set URL_PREFIX=ftp://ftp.ncbi.nlm.nih.gov.
    os.environ["URL_PREFIX"] = "ftp://ftp.ncbi.nlm.nih.gov"

    ## Instances with adequate resources for making gsnap and rapsearch indexes.
    # Need to have write access to destination folders.
    os.environ["RAPSEARCH_SERVER_IP"] = "54.191.193.210"
    os.environ["GSNAP_SERVER_IP"] = "34.211.67.166"

    ##### COMMANDS #####
    date_setup()
    make_gsnap_index()
    make_rapsearch_index()
    make_lineage_files()
    make_accession_mapping()


def make_gsnap_index():
    nt = "/blast/db/FASTA/nt.gz"
    date = set_index_date([nt])
    dest = "s3://idseq-database/alignment_indexes/" + date
    os.environ["INPUT_FASTA_S3"] = os.environ["URL_PREFIX"] + nt
    os.environ["SERVER_IP"] = os.environ["GSNAP_SERVER_IP"]
    os.environ["KEY_S3_PATH"] = "s3://idseq-secrets/idseq-production.pem"
    os.environ["OUTPUT_PATH_S3"] = dest
    os.environ["OUTPUT_NAME"] = "nt_k16"
    execute_command("idseq_pipeline gsnap_indexing")


def make_rapsearch_index():
    nr = "/blast/db/FASTA/nr.gz"
    date = set_index_date([nr])
    dest = "s3://idseq-database/alignment_indexes/" + date
    os.environ["INPUT_FASTA_S3"] = os.environ["URL_PREFIX"] + nr
    os.environ["SERVER_IP"] = os.environ["RAPSEARCH_SERVER_IP"]
    os.environ["KEY_S3_PATH"] = "s3://czbiohub-infectious-disease/idseq-alpha.pem"
    os.environ["OUTPUT_PATH_S3"] = dest
    os.environ["OUTPUT_NAME"] = "nr_rapsearch"
    pipeline_command("rapsearch_indexing")


def make_lineage_files():
    taxdump = "/pub/taxonomy/taxdump.tar.gz"
    date = set_index_date([taxdump])
    dest = "s3://idseq-database/taxonomy/" + date
    os.environ["OUTPUT_PATH_S3"] = dest
    os.environ["INPUT"] = os.environ["URL_PREFIX"] + taxdump
    pipeline_command("lineages")


def make_accession_mapping():
    nr = "/blast/db/FASTA/nr.gz"
    nt = "/blast/db/FASTA/nt.gz"
    mapping_files = ["/pub/taxonomy/accession2taxid/nucl_est.accession2taxid.gz",
                     "/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
                     "/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz",
                     "/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
                     "/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz",
                     "/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"]
    input_files = [nr] + [nt] + mapping_files
    date = set_index_date(input_files)
    pre = os.environ["URL_PREFIX"]
    dest = "s3://idseq-database/alignment_data/" + date
    map_env = [pre + f for f in mapping_files]
    map_env = ','.join(map_env)  # Concatenate into env var

    os.environ["MAPPING_FILES"] = map_env
    prev_mapping = "s3://czbiohub-infectious-disease/references/accession2taxid.db.gz"
    cmd = "curate_accession2taxid --mapping_files %s --nr_file %s --nt_file %s --output_s3_folder %s --previous_mapping %s".format(
        map_env, pre + nr, pre + nt, dest, prev_mapping)
    pipeline_command(cmd)


def date_setup():
    execute_command("mkdir -p %s" % local_work_dir)
    execute_command("aws s3 cp --quiet %s %s/" % (ncbitool_s3_path, local_work_dir))
    execute_command("chmod u+x %s/ncbitool" % local_work_dir)


def set_index_date(input_files):
    ncbitool_path = local_work_dir + "/ncbitool"
    maxe = None
    for f in input_files:
        # Time from ncbitool is already in UTC. Get max date in input file set.
        command = "%s file %s" % (ncbitool_path, f)
        output = execute_command_with_output(command).split("File Info: ")[1]
        date = json.loads(output)["ModTime"]
        date = datetime.datetime.strptime(date, "%Y-%m-%dT%H:%M:%S")
        if maxe is None or date > maxe:
            maxe = date

    # Format source file timestamp
    source_time = datetime.datetime.strftime(maxe, "%Y-%m-%d-utc-")
    utime = str(int(time.mktime(maxe.timetuple())))
    source_str = source_time + utime + "-unixtime"

    # Format current timestamp
    now = datetime.datetime.utcnow()
    utime = str(int(time.time()))
    now_str = datetime.datetime.strftime(now, "%Y-%m-%d-utc-") + utime + "-unixtime"

    res = source_str + "__" + now_str
    return res


def pipeline_command(cmd):
    cmd = cmd.split()
    proc = subprocess.Popen(["idseq_pipeline"] + cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    res = proc.communicate()
    for ln in (res[0] + res[1]).splitlines():
        # Print out stdout and stderr
        print ln


main()
