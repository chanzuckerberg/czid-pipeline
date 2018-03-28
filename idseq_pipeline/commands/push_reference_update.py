"""Command for pushing a reference update."""
import json
import os

from .base import Base
from .common import env_set_if_blank, execute_command, install_ncbitool_locally, execute_command_with_output, \
    datetime as dt, time

LOCAL_WORK_DIR = "idseq_pipeline_temp"


class Push_reference_update(Base):
    def run(self):
        ##### PARAMETERS TO EDIT #####

        ## To pull NCBI references through ncbitool, which archives and records versions, set URL_PREFIX=''.
        # To pull the latest files from NCBI without using ncbitool (no version information), set URL_PREFIX=ftp://ftp.ncbi.nlm.nih.gov.
        env_set_if_blank("URL_PREFIX", "ftp://ftp.ncbi.nlm.nih.gov")

        ## Instances with adequate resources for making gsnap and rapsearch indexes.
        # Need to have write access to destination folders.
        env_set_if_blank("RAPSEARCH_SERVER_IP", "54.191.193.210")
        env_set_if_blank("GSNAP_SERVER_IP", "34.211.67.166")
        env_set_if_blank("DEST_PREFIX", "s3://idseq-database")

        ##### COMMANDS #####
        install_ncbitool_locally(LOCAL_WORK_DIR)
        self.make_gsnap_index()
        self.make_rapsearch_index()
        self.make_accession_mapping()

    def make_gsnap_index(self):
        nt = "/blast/db/FASTA/nt.gz"
        date = self.set_index_date([nt])
        dest = os.environ["DEST_PREFIX"] + "/alignment_indexes/" + date
        os.environ["INPUT_FASTA_S3"] = os.environ["URL_PREFIX"] + nt
        os.environ["SERVER_IP"] = os.environ["GSNAP_SERVER_IP"]
        os.environ["KEY_S3_PATH"] = "s3://idseq-secrets/idseq-production.pem"
        os.environ["OUTPUT_PATH_S3"] = dest
        os.environ["OUTPUT_NAME"] = "nt_k16"
        execute_command("idseq_pipeline gsnap_indexing")

    def make_rapsearch_index(self):
        nr = "/blast/db/FASTA/nr.gz"
        date = self.set_index_date([nr])
        dest = os.environ["DEST_PREFIX"] + "/alignment_indexes/" + date
        os.environ["INPUT_FASTA_S3"] = os.environ["URL_PREFIX"] + nr
        os.environ["SERVER_IP"] = os.environ["RAPSEARCH_SERVER_IP"]
        os.environ["KEY_S3_PATH"] = "s3://czbiohub-infectious-disease/idseq-alpha.pem"
        os.environ["OUTPUT_PATH_S3"] = dest
        os.environ["OUTPUT_NAME"] = "nr_rapsearch"
        execute_command("idseq_pipeline rapsearch_indexing")

    def make_accession_mapping(self):
        nr = "/blast/db/FASTA/nr.gz"
        nt = "/blast/db/FASTA/nt.gz"
        mapping_files = ["/pub/taxonomy/accession2taxid/nucl_est.accession2taxid.gz",
                         "/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
                         "/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz",
                         "/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
                         "/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz",
                         "/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"]
        input_files = [nr] + [nt] + mapping_files
        date = self.set_index_date(input_files)
        pre = os.environ["URL_PREFIX"]
        dest = os.environ["DEST_PREFIX"] + "/alignment_data/" + date
        map_env = [pre + f for f in mapping_files]
        map_env = ','.join(map_env)  # Concatenate into env var

        os.environ["MAPPING_FILES"] = map_env
        prev_mapping = "s3://czbiohub-infectious-disease/references/accession2taxid.db.gz"
        cmd = "idseq_pipeline curate_accession2taxid --mapping_files %s --nr_file %s --nt_file %s --output_s3_folder %s --previous_mapping %s".format(
            map_env, pre + nr, pre + nt, dest, prev_mapping)
        execute_command(cmd)

    def set_index_date(self, input_files):
        ncbitool_path = LOCAL_WORK_DIR + "/ncbitool"
        maxe = None
        for f in input_files:
            # Time from ncbitool is already in UTC. Get max date in input file set.
            command = "%s file %s" % (ncbitool_path, f)
            output = execute_command_with_output(command).split("File Info: ")[1]
            date = json.loads(output)["ModTime"]
            date = dt.datetime.strptime(date, "%Y-%m-%dT%H:%M:%S")
            if maxe is None or date > maxe:
                maxe = date

        # Format source file timestamp
        source_time = dt.datetime.strftime(maxe, "%Y-%m-%d-utc-")
        utime = str(int(time.mktime(maxe.timetuple())))
        source_str = source_time + utime + "-unixtime"

        # Format current timestamp
        now = dt.datetime.utcnow()
        utime = str(int(time.time()))
        now_str = dt.datetime.strftime(now, "%Y-%m-%d-utc-") + utime + "-unixtime"

        res = source_str + "__" + now_str
        return res
