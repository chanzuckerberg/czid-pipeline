import shelve
import re
import json
import threading
import traceback
import os

import subprocess
try:
    from subprocess import DEVNULL #pylint: disable=no-name-in-module
except:
    DEVNULL = open(os.devnull, "r+b")

from .common import * #pylint: disable=wildcard-import

REF_DISPLAY_RANGE = 100
MAX_SEQ_DISPLAY_SIZE = 6000


# Test with the following function call
# generate_alignment_viz_json('../../nt','nt.db','NT', 'taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8', 'taxid_annot_sorted_nt.fasta', 'align_viz')

def parse_reads(annotated_fasta, db_type):
    read2seq = {}

    search_string = "species_%s" % (db_type.lower())
    adv_search_string = "family_%s:([-\d]+):.*genus_%s:([-\d]+):.*species_%s:([-\d]+).*NT:[^:]*:(.*)" % (db_type.lower(), db_type.lower(), db_type.lower())

    with open(annotated_fasta, 'r') as af:
        read_id = ''
        sequence = ''
        for line in af:
            if line[0] == '>':
                read_id = line
            else:
                sequence = line
                m = re.search("%s:([\d-]*)" % search_string, read_id)
                if m:
                    if int(m.group(1)) > 0: # it's a match
                        ma = re.search(adv_search_string, read_id)
                        if ma:
                            read2seq[ma.group(4).rstrip()] = [sequence.rstrip(), ma.group(1), ma.group(2), ma.group(3)]
    return read2seq


def generate_alignment_viz_json(nt_file, nt_loc_db, db_type,
                                annotated_m8, annotated_fasta,
                                output_json_dir):
    """Generate alignment details from the reference sequence, m8 and annotated fasta """
    # Go through annotated_fasta with a db_type (NT/NR match).  Infer the family/genus/species info

    if db_type != 'NT' and db_type != 'NR':
        return

    read2seq = parse_reads(annotated_fasta, db_type)

    print("Read to Seq dictionary size: %d" % len(read2seq))

    # Go through m8 file infer the alignment info. grab the fasta sequence, lineage info
    groups = {}
    line_count = 0
    nt_loc_dict = shelve.open(nt_loc_db)
    with open(annotated_m8, 'r') as m8f:
        for line in m8f:
            line_count += 1
            if line_count % 100000 == 0:
                print("%d lines in the m8 file processed." % line_count)
            line_columns = line.rstrip().split("\t")
            read_id = extract_m8_readid(line_columns[0])
            seq_info = read2seq.get(read_id)
            if seq_info:
                accession_id = line_columns[1]
                metrics = line_columns[2:]
                # "ad" is short for "accession_dict" aka "accession_info"
                ad = groups.get(accession_id, {'reads': []})
                sequence, ad['family_id'], ad['genus_id'], ad['species_id'] = seq_info
                ref_start = int(metrics[-4])
                ref_end = int(metrics[-3])
                if ref_start > ref_end: # SWAP
                    (ref_start, ref_end) = (ref_end, ref_start)
                ref_start -= 1
                prev_start = (ref_start - REF_DISPLAY_RANGE) if (ref_start - REF_DISPLAY_RANGE) > 0 else 0
                post_end = ref_end + REF_DISPLAY_RANGE
                ad['reads'].append([read_id, sequence, metrics, (prev_start, ref_start, ref_end, post_end)])
                ad['ref_link'] = "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=fasta" % accession_id
                groups[accession_id] = ad

    print("%d lines in the m8 file" % line_count)
    print("%d unique accession ids" % len(groups))

    if nt_file.startswith("s3://"):
        get_sequences_by_accession_list_from_s3(groups, nt_loc_dict, nt_file)
    else:
        get_sequences_by_accession_list_from_file(groups, nt_loc_dict, nt_file)

    result_dict = {}
    to_be_deleted = []
    for accession_id, ad in groups.iteritems():
        tmp_file = 'accession-%s' % accession_id
        if ad['ref_seq_len'] <= MAX_SEQ_DISPLAY_SIZE and 'ref_seq' not in ad:
            with open(tmp_file, "rb") as tf:
                ad['ref_seq'] = tf.read()
                # the previous value of ad['ref_seq_len'] is greater because it includes newline characters
                # the value we set here excludes newline characters
                ad['ref_seq_len'] = len(ad['ref_seq'])
            to_be_deleted.append(tmp_file)
        if 'ref_seq' in ad:
            ref_seq = ad['ref_seq']
            for read in ad['reads']:
                prev_start, ref_start, ref_end, post_end = read[3]
                read[3] = [ref_seq[prev_start:ref_start], ref_seq[ref_start:ref_end], ref_seq[ref_end:post_end]]
        else:
            # the reference sequence is too long to read entirely in RAM, so we only read the mapped segments
            tmp_file = 'accession-%s' % accession_id
            with open(tmp_file, "rb") as tf:
                for read in ad['reads']:
                    prev_start, ref_start, ref_end, post_end = read[3]
                    tf.seek(prev_start, 0)
                    segment = tf.read(post_end - prev_start)
                    read[3] = [
                        segment[0:(ref_start - prev_start)],
                        segment[(ref_start - prev_start):(ref_end - prev_start)],
                        segment[(ref_end - prev_start):(post_end - prev_start)]
                    ]
            to_be_deleted.append(tmp_file)
        if ad['ref_seq_len'] > MAX_SEQ_DISPLAY_SIZE:
            ad['ref_seq'] = '...Reference Seq Too Long ...'
        family_id = ad.pop('family_id')
        genus_id = ad.pop('genus_id')
        species_id = ad.pop('species_id')
        family_dict = result_dict.get(family_id, {})
        genus_dict = family_dict.get(genus_id, {})
        species_dict = genus_dict.get(species_id, {})
        species_dict[accession_id] = ad
        genus_dict[species_id] = species_dict
        family_dict[genus_id] = genus_dict
        result_dict[family_id] = family_dict

    deleter_thread = threading.Thread(target=map, args=[os.remove, to_be_deleted])
    deleter_thread.start()

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

    deleter_thread.join()

    summary = "Read2Seq Size: %d, M8 lines %d, %d unique accession ids" % (len(read2seq), line_count, len(groups))
    summary_file_name = "%s.summary" % output_json_dir
    with open(summary_file_name, 'w') as summaryf:
        summaryf.write(summary)
    return summary_file_name


def delete_many(files, semaphore=None): #pylint: disable=dangerous-default-value
    try:
        for f in files:
            os.remove(f)
    except:
        with print_lock:
            print("Couldn't delete some temp files.  Moving on.")
    finally:
        if semaphore:
            semaphore.release()

def get_sequences_by_accession_list_from_file(accession2seq, nt_loc_dict, nt_file):
    with open(nt_file) as ntf:
        for accession_id, accession_info in accession2seq.iteritems():
            accession_info['ref_seq'] = get_sequence_by_accession_id_ntf(accession_id, nt_loc_dict, ntf)
            accession_info['ref_seq_len'] = len(accession_info['ref_seq'])

def get_sequences_by_accession_list_from_s3(accession_id_groups, nt_loc_dict, nt_s3_path):
    threads = []
    error_flags = {}
    semaphore = threading.Semaphore(64)
    mutex = threading.RLock()
    nt_bucket, nt_key = nt_s3_path[5:].split("/", 1)
    for accession_id, accession_info in accession_id_groups.iteritems():
        semaphore.acquire()
        t = threading.Thread(
            target=get_sequence_for_thread,
            args=[error_flags, accession_info, accession_id, nt_loc_dict, nt_bucket, nt_key, semaphore, mutex]
        )
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    if error_flags:
        raise "Sorry there was an error"

def get_sequence_for_thread(error_flags, accession_info, accession_id, nt_loc_dict, nt_bucket, nt_key, semaphore, mutex, seq_count=[0]): #pylint: disable=dangerous-default-value
    try:
        ref_seq_len = get_sequence_by_accession_id_s3(accession_id, nt_loc_dict, nt_bucket, nt_key)
        with mutex:
            accession_info['ref_seq_len'] = ref_seq_len
            seq_count[0] += 1
            if seq_count[0] % 100 == 0:
                print("%d sequences fetched, most recently %s" % (seq_count[0], accession_id))
    except:
        with mutex:
            if not error_flags:
                traceback.print_exc()
            error_flags["error"] = 1
    finally:
        semaphore.release()

def get_sequence_by_accession_id_ntf(accession_id, nt_loc_dict, ntf):
    ref_seq = ''
    entry = nt_loc_dict.get(accession_id)
    if entry:
        range_start = entry[0]+entry[1]
        seq_len = entry[2]
        ntf.seek(range_start, 0)
        ref_seq = ntf.read(seq_len).replace("\n", "")
    return ref_seq

def get_sequence_by_accession_id_s3(accession_id, nt_loc_dict, nt_bucket, nt_key):
    seq_len = 0
    entry = nt_loc_dict.get(accession_id)
    if entry:
        range_start = entry[0]+entry[1]
        seq_len = entry[2]
        pipe_file = 'pipe-accession-%s' % accession_id
        accession_file = 'accession-%s' % accession_id
        os.mkfifo(pipe_file)
        get_range = "aws s3api get-object --range bytes=%d-%d --bucket %s --key %s %s" % (range_start, range_start + seq_len - 1, nt_bucket, nt_key, pipe_file)
        get_range_proc = subprocess.Popen(get_range, shell=True, stdout=DEVNULL)
        subprocess.check_call("tr -d '\n' < {pipe_file} > {accession_file}".format(pipe_file=pipe_file, accession_file=accession_file), shell=True)
        exitcode = get_range_proc.wait()
        assert exitcode == 0
        os.remove(pipe_file)
    return seq_len

def accessionid2seq_main(arguments):

    # Make work directory
    dest_dir = os.path.join(DEST_DIR, "accession2seq/tmp-%d" % os.getpid())
    execute_command("mkdir -p %s" % dest_dir)

    s3_db_path = arguments.get('--s3_db_path')
    s3_db_loc_path = arguments.get('--s3_db_loc_path')
    db_type = arguments.get('--db_type')
    input_fasta_s3_path = arguments.get('--input_fasta_s3_path')
    input_m8_s3_path = arguments.get('--input_m8_s3_path')
    output_json_s3_path = arguments.get('--output_json_s3_path')
    output_json_s3_path = arguments.get('--output_json_s3_path').rstrip('/')

    db_path = arguments.get('--local_db_path') # Try to get local file first


    local_db_loc_path = os.path.join(dest_dir, os.path.basename(s3_db_loc_path))
    local_fasta_path = os.path.join(dest_dir, os.path.basename(input_fasta_s3_path))
    local_m8_path = os.path.join(dest_dir, os.path.basename(input_m8_s3_path))
    local_json_path = os.path.join(dest_dir, "align_viz")

    if not db_path:
        db_path = s3_db_path
    execute_command("aws s3 cp --quiet %s %s" %(s3_db_loc_path, local_db_loc_path))
    execute_command("aws s3 cp --quiet %s %s" %(input_fasta_s3_path, local_fasta_path))
    execute_command("aws s3 cp --quiet %s %s" %(input_m8_s3_path, local_m8_path))
    summary_file_name = generate_alignment_viz_json(db_path, local_db_loc_path, db_type,
                                                    local_m8_path, local_fasta_path, local_json_path)
    # copy the data over
    execute_command("aws s3 cp --quiet %s %s --recursive" % (local_json_path, output_json_s3_path))
    execute_command("aws s3 cp --quiet %s %s/" % (summary_file_name, os.path.dirname(output_json_s3_path)))

    # Clean up
    execute_command("rm -rf %s" % dest_dir)
