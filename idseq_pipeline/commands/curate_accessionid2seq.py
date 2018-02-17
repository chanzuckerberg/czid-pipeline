"""The curate_accessionid2seq command."""
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

def curate_seq_dict(nt_file, output_db_file):
    """Curate nt/nr accession id to  file offset to get the sequence"""
    print "Read the nt/nr file"
    # Read the nt/nr file
    taxon_map = shelve.open(output_db_file)
    with open(nt_file) as ntf:
        seq_offset = 0
        seq_len = 0
        header_len = 0
        lines = 0
        accession_id = ""
        for line in ntf:
            lines += 1
            if lines % 100000 == 0:
                print "%3.1f M lines. " %  (lines/1000000.0)
            if line[0] == '>': # header line
                if seq_len > 0 and len(accession_id) > 0:
                    taxon_map[accession_id] = (seq_offset, header_len, seq_len)
                seq_offset = seq_offset + header_len + seq_len
                header_len = len(line)
                seq_len = 0
                s = re.match('^>([^ ]*).*', line)
                if s:
                    accession_id = s.group(1)
            else:
                seq_len += len(line)
        if seq_len > 0 and len(accession_id) > 0:
            taxon_map[accession_id] = (seq_offset, header_len, seq_len)
    taxon_map.close()

class Curate_accessionid2seq(Base):
    def run(self):
        from .common import * #TO DO: clean up the imports across this package
