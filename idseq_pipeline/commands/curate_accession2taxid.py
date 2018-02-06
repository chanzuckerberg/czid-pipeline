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

def curate_taxon_dict(results):
    print "Read the nt/nr file"
    # Read the nt/nr file
    dbids = {}
    lines = 0
    for gzf in [results['--nr_file'], results['--nt_file']]:  
        with gzip.open(gzf, 'rb') as seqf:
            for line in seqf:
                lines += 1
                if lines % 100000 == 0:
                    print "%f M lines. " %  (lines/1000000.0)
                if line[0] == '>': # header line
                    s = re.match('^>([^ \.]*).*', line)
                    if s:
                        dbids[s.group(1)] = 1
    print "Read the mapping file and select ..."
    # Read the accession2tax mappng file
    outf = open(results['--output_mapping_file'], 'wb')
    lines = 0
    with gzip.open(results['--mapping_file'], 'rb') as mapf:
        for line in mapf:
            fields = line.split("\t")
            lines += 1
            if lines % 100000 == 0:
                print "%f M lines. %s, %s" %  (lines/1000000.0, fields[0], fields[2])
            if dbids.get(fields[0]):
                outf.write(line)
    outf.close()

class Curate_accession2taxid(Base):
    """Curate dictionary based on existence in NT/NR and accessiont2taxon mapping"""

    def run(self):
        #to do: get the file from ncbitool
        curate_taxon_dict(self.options)
        #to do: generate_db()
