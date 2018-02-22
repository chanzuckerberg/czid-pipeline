"""The  accessionid2seq command."""
from .base import Base
from .common import * #pylint: disable=wildcard-import
import sys
import time
import subprocess
import random
import pickle
import shelve
import argparse
import gzip
import re
import json
import threading
import traceback


# Test with the following function call
# generate_alignment_viz_json('../../nt','nt.db','NT', 'taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8', 'taxid_annot_sorted_nt.fasta', 'align_viz')

def generate_alignment_viz_json(nt_file, nt_loc_db, db_type,
                                annotated_m8, annotated_fasta,
                                output_json_dir):
    # Go through annotated_fasta with a db_type (NT/NR match). Infer the family/genus/species info
    if db_type != 'NT' and db_type != 'NR':
        return

    read2seq = {}

    search_string = "species_%s" % (db_type.lower())
    adv_search_string ="family_%s:([-\d]+):.*genus_%s:([-\d]+):.*species_%s:([-\d]+).*NT:[^:]*:(.*)" % (db_type.lower(), db_type.lower(), db_type.lower())

    with open(annotated_fasta, 'r') as af:
        read_id = ''
        sequence = ''
        for line in af:
            if line[0] == '>':
                read_id = line
            else:
                sequence = line
                m = re.search("%s:([\d-]*)" % search_string,read_id)
                if m:
                    if int(m.group(1)) > 0: # it's a match
                        ma = re.search(adv_search_string, read_id)
                        if ma:
                            read2seq[ma.group(4).rstrip()] = [sequence.rstrip(), ma.group(1), ma.group(2), ma.group(3)]
    print("Read to Seq dictionary size: %d" % len(read2seq))

    # Go through m8 file infer the alignment info. grab the fasta sequence, lineage info, construct a tree, get the actual sequence from nt_file
    nt_loc_dic = shelve.open(nt_loc_db)
    nt = open(nt_file)
    result_dict = {}
    line_count = 0
    with open(annotated_m8, 'r') as m8f:
        for line in m8f:
            line_count += 1
            line_columns = line.rstrip().split("\t")
            read_id = extract_m8_readid(line_columns[0])
            accession_id = line_columns[1]
            metrics = line_columns[2:]
            seq_info = read2seq.get(read_id)
            if seq_info:
                [sequence, family_id, genus_id, species_id] = seq_info
                result_dict[family_id] = result_dict.get(family_id, {})
                result_dict[family_id][genus_id] = result_dict[family_id].get(genus_id, {})
                result_dict[family_id][genus_id][species_id] = result_dict[family_id][genus_id].get(species_id, {})
                accession_dict = result_dict[family_id][genus_id][species_id].get(accession_id, {})
                accession_dict['reads'] = accession_dict.get('reads', [])
                accession_dict['reads'].append((read_id, sequence, metrics))
                if not accession_dict.get('ref_seq'):
                    # get reference sequence
                    entry = nt_loc_dic.get(accession_id)
                    if entry:
                        nt.seek(entry[0]+entry[1], 0)
                        ref_seq = nt.read(entry[2])
                        accession_dict['ref_seq'] = ref_seq.replace("\n", "")
                    else:
                        accession_dict['ref_seq'] = 'NOT FOUND'

                result_dict[family_id][genus_id][species_id][accession_id] = accession_dict
    print("%d lines in the m8 file" % line_count")

    # output json by species, genus, family
    execute_command("mkdir -p %s " % output_json_dir)
    for (family_id, family_dict) in result_dict.iteritems():
        with open("%s/%s.family.%d.align_viz.json" %(output_json_dir, db_type.lower(), int(family_id)), 'wb') as outjf:
            json.dump(family_dict, outjf)
        for (genus_id, genus_dict) in family_dict.iteritems():
            with open("%s/%s.genus.%d.align_viz.json" %(output_json_dir, db_type.lower(), int(genus_id)), 'wb') as outjf:
                json.dump(genus_dict, outjf)
            for (species_id, species_dict) in genus_dict.iteritems():
                with open("%s/%s.species.%d.align_viz.json" %(output_json_dir, db_type.lower(), int(species_id)), 'wb') as outjf:
                    json.dump(species_dict, outjf)

class Accessionid2seq(Base):
    def run(self):
        print("Placeholder")
