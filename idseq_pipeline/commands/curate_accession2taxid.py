"""The curate_accession2taxid command."""
from .base import Base
import sys
import time
import subprocess
import random
import pickle
import shelve
import argparse
import gzip
import re
import threading
import traceback

def curate_taxon_dict(nt_file, nr_file, mapping_files, output_mapping_file):
    """Curate accessiont2taxid mapping based on existence in NT/NR"""
    print "Read the nt/nr file"
    # Read the nt/nr file
    dbids = set()
    lines = 0
    for gzf in [nr_file, nt_file]:
        with gzip.open(gzf, 'rb') as seqf:
            for line in seqf:
                lines += 1
                if lines % 100000 == 0:
                    print "%3.1f M lines. " %  (lines/1000000.0)
                if line[0] == '>': # header line
                    s = re.match('^>([^ \.]*).*', line)
                    if s:
                        dbids.add(s.group(1))
    print "Read the mapping file and select ..."
    # Read the accession2taxid mapping files
    outf = open(output_mapping_file, 'wb')
    lines = 0
    for mapping_file in mapping_files:
        with gzip.open(mapping_file, 'rb') as mapf:
            for line in mapf:
                fields = line.split("\t")
                lines += 1
                if lines % 100000 == 0:
                    print "%3.1f M lines. %s, %s" %  (lines/1000000.0, fields[0], fields[2])
                if fields[0] in dbids:
                    outf.write(line)
    outf.close()

def generate_accession2taxid_db(mapping_file, output_db_file, input_gzipped):
    lines = 0
    taxon_map = shelve.open(output_db_file)
    if input_gzipped:
        mapf = gzip.open(mapping_file, 'rb')
    else:
        mapf = open(mapping_file, 'rb')
    for line in mapf:
        fields = line.rstrip().split("\t")
        if len(fields) == 4:
            lines += 1
            if lines % 100000 == 0:
                print "%d lines. %s, %s" %  (lines, fields[0], fields[2])
            accession_id = fields[0]
            taxon_id = fields[2]
            taxon_map[accession_id] = taxon_id
    print "close the db file"
    taxon_map.close()
    mapf.close()

class Curate_accession2taxid(Base):
    def run(self):
        from .common import * #TO DO: clean up the imports across this package

        # Make work directory
        dest_dir = os.path.join(DEST_DIR, 'accession2taxid')
        execute_command("mkdir -p %s" % dest_dir)

        # Retrieve the reference files
        print "Retrieving references"
        arguments = self.options
        threads = {}
        threads['nt'] = MyThread(target=download_reference_locally_with_version_any_source_type,
            args=(arguments['--nt_file'], dest_dir, dest_dir))
        threads['nt'].start()
        threads['nr'] = MyThread(download_reference_locally_with_version_any_source_type,
            (arguments['--nr_file'], dest_dir, dest_dir))
        threads['nr'].start()
        mapping_files_local = []
        mapping_version_numbers = []
        mapping_files_sources = arguments['--mapping_files'].split(",")
        mapping_file_results = {}
        for f in mapping_files_sources:
            assert f not in threads
            threads[f] = MyThread(target=download_reference_locally_with_version_any_source_type,
                args=(f, dest_dir, dest_dir))
            threads[f].start()
        for f in threads:
            threads[f].join()
            assert not threads[f].exception
        nt_file_local, nt_version_number = threads['nt'].result
        nr_file_local, nr_version_number = threads['nr'].result
        for f in mapping_files_sources:
            mapping_file_local, mapping_version_number = threads[f].result
            mapping_files_local.append(mapping_file_local)
            mapping_version_numbers.append(mapping_version_number)
        print "Reference download finished"

        # Produce the output file
        print "Curating accession2taxid"
        output_mapping_file = os.path.join(DEST_DIR, 'curated_accession2taxid.txt')
        curate_taxon_dict(nt_file_local, nr_file_local, mapping_files_local, output_mapping_file)
        print "Curation finished"

        # Convert to a berkeley db
        print "Writing berkeley db"
        output_db_file = os.path.join(dest_dir, 'accession2taxid.db')
        execute_command("rm -f %s" % output_db_file)
        generate_accession2taxid_db(output_mapping_file, output_db_file, False)
        output_s3_path = arguments['--output_s3_folder'].rstrip('/')
        output_s3_file_gz = os.path.join(output_s3_path, os.path.basename(output_db_file) + ".gz")

        # Record difference with old accession2taxid
        previous_mapping_s3 = arguments.get('--previous_mapping')
        if previous_mapping_s3:
            print "Computing diff"
            previous_mapping_local = os.path.join(dest_dir, "previous_mapping.db.gz")
            execute_command("aws s3 cp --quiet %s %s" % (previous_mapping_s3, previous_mapping_local))
            print "Starting gunzip"
            execute_command("gunzip -fk %s" % previous_mapping_local)
            print "Starting shelve open"
            previous_mapping = shelve.open(os.path.splitext(previous_mapping_local)[0])
            new_mapping = shelve.open(output_db_file)
            added_accessionids = []
            removed_accessionids = []
            print "Starting to diff mappings"
            for new_key in new_mapping:
                if new_key not in previous_mapping:
                    added_accessionids.append(new_key)
            for previous_key in previous_mapping:
                if previous_key not in new_mapping:
                    removed_accessionids.append(previous_key)
            diff_accessionids = { 'old_file': previous_mapping_s3, 'new_file': output_s3_file_gz,
                                  'added': added_accessionids, 'removed': removed_accessionids }
            diff_accessionids_file = os.path.join(dest_dir, "accession_diff.txt")
            print "Starting json dump"
            with open(diff_accessionids_file, 'wb') as f:
                json.dump(diff_accessionids, f)
            execute_command("aws s3 cp --quiet %s %s/" % (diff_accessionids_file, output_s3_path))

        # Upload result and record versions
        print "Uploading result to S3"
        execute_command("gzip -c {output_db_file} | aws s3 cp --quiet - {output_s3_file_gz}".format(output_db_file=output_db_file, output_s3_file_gz=output_s3_file_gz))
        upload_version_tracker(mapping_files_sources + [arguments['--nt_file'], arguments['--nr_file']],
                               'accession2taxid',
                               mapping_version_numbers + [nt_version_number, nr_version_number],
                               output_s3_path, self.version)
