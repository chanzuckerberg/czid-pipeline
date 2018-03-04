import shelve
import re
import json
import threading
import multiprocessing
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

def get_chunk_start(full_size, num_chunks, chunk_idx):
    if chunk_idx >= num_chunks:
        return full_size
    return int(chunk_idx * (full_size / float(num_chunks)))

def get_dict_chunk(entire_dict, num_chunks, chunk_idx):
    chunk_start = get_chunk_start(len(entire_dict), num_chunks, chunk_idx)
    chunk_end = get_chunk_start(len(entire_dict), num_chunks, chunk_idx + 1)
    result = {}
    all_keys = sorted(entire_dict.keys())
    for k in all_keys[chunk_start:chunk_end]:
        result[k] = entire_dict[k]
    return result

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
                                output_json_dir,
                                slave_id=None,
                                num_slaves=8):
    """Generate alignment details from the reference sequence, m8 and annotated fasta """
    # Go through annotated_fasta with a db_type (NT/NR match).  Infer the family/genus/species info

    if db_type != 'NT' and db_type != 'NR':
        return

    if slave_id == None:  # master
        tentative_summary = "%s.summary.tentative" % output_json_dir
        final_summary = "%s.summary" % output_json_dir
        if os.path.exists(final_summary):
            os.remove(final_summary)
        if os.path.exists(tentative_summary):
            os.remove(tentative_summary)
        execute_command("mkdir -p %s " % output_json_dir)
        processes = [multiprocessing.Process(target=generate_alignment_viz_json,
                                             args=[nt_file, nt_loc_db, db_type, annotated_m8, annotated_fasta,
                                                   output_json_dir, sid, num_slaves])
                     for sid in range(num_slaves)]
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        for p in processes:
            if p.exitcode != 0:
                raise "Sorry, a generate_alignment_viz_json subprocess encountered a problem."
        os.rename(tentative_summary, final_summary)
        return final_summary

    read2seq = parse_reads(annotated_fasta, db_type)

    if slave_id == 0:
        print("Read to Seq dictionary size: %d" % len(read2seq))

    # Go through m8 file infer the alignment info. grab the fasta sequence, lineage info, construct a tree, get the actual sequence from nt_file
    nt_loc_dict = shelve.open(nt_loc_db)
    accession2seq = {}
    with open(annotated_m8, 'r') as m8f:
        for line in m8f:
            line_columns = line.rstrip().split("\t")
            accession_id = line_columns[1]
            accession2seq[accession_id] = accession2seq.get(accession_id, {'count': 0})
            accession2seq[accession_id]['count'] += 1
    if slave_id == 0:
        print("%d unique accession ids" % len(accession2seq))

    # Replace accession2seq with the part this slave will work on.
    # TODO partition accession ids based on range lengths so we balance across slaves
    accession_ids_count = len(accession2seq)
    accession2seq = get_dict_chunk(accession2seq, num_slaves, slave_id)

    if nt_file.startswith("s3://"):
        # This can take an hour easily.
        get_sequences_by_accession_list_from_s3(accession2seq, nt_loc_dict, nt_file, num_slaves, slave_id)
    else:
        get_sequences_by_accession_list_from_file(accession2seq, nt_loc_dict, nt_file)

    # group reads by accession id
    slave_results = {}
    line_count = 0

    with open(annotated_m8, 'r') as m8f:
        for line in m8f:
            line_count += 1
            line_columns = line.rstrip().split("\t")
            accession_id = line_columns[1]
            if accession_id in accession2seq:
                ref_seq = accession2seq[accession_id]['seq']
                read_id = extract_m8_readid(line_columns[0])
                seq_info = read2seq.get(read_id)
                metrics = line_columns[2:]
                if slave_id == 0 and line_count % 100000 == 0:
                    print("%d lines in the m8 file processed by slave 0." % line_count)

                if seq_info:
                    if accession_id not in slave_results:
                        slave_results[accession_id] = {'accession_id': accession_id, 'reads': []}
                    ad = slave_results[accession_id]
                    sequence, ad['family_id'], ad['genus_id'], ad['species_id'] = seq_info
                    ref_start = int(metrics[-4])
                    ref_end = int(metrics[-3])
                    if ref_start > ref_end: # SWAP
                        (ref_start, ref_end) = (ref_end, ref_start)
                    ref_start -= 1
                    prev_start = (ref_start - REF_DISPLAY_RANGE) if (ref_start - REF_DISPLAY_RANGE) > 0 else 0
                    post_end = ref_end + REF_DISPLAY_RANGE
                    if len(ref_seq) > MAX_SEQ_DISPLAY_SIZE:
                        ad['ref_seq'] = '...Reference Seq Too Long ...'
                    else:
                        ad['ref_seq'] = ref_seq
                    ref_read = [ref_seq[prev_start:ref_start], ref_seq[ref_start:ref_end], ref_seq[ref_end:post_end]]
                    ad['reads'].append((read_id, sequence, metrics, ref_read))
                    ad['ref_seq_len'] = len(ref_seq)
                    ad['ref_link'] = "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=fasta" % accession_id

    if slave_id == 0:
        print("%d lines in the m8 file" % line_count)

    for accession_id, ad in slave_results:
        with open("%s/%s__family__%d__genus__%d__species__%d__accession__%s__align_viz.json" %(output_json_dir, db_type.lower(), int(ad['family_id']), int(ad['genus_id']), int(ad['species_id']), accession_id), 'wb') as outjf:
            json.dump(ad, outjf)

    if slave_id == 0:
        summary = "Read2Seq Size: %d, M8 lines %d, %d unique accession ids" % (len(read2seq), line_count, accession_ids_count)
        summary_file_name = "%s.summary.tentative" % output_json_dir
        with open(summary_file_name, 'w') as summaryf:
            summaryf.write(summary)



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

def async_delete(tmp_file, mutex=threading.RLock(), to_be_deleted=[[]], batch_size=1024, high_watermark=2048, max_threads=threading.Semaphore(3)): #pylint: disable=dangerous-default-value
    with mutex:
        tbd = to_be_deleted[0]
        if tmp_file:
            tbd.append(tmp_file)
            watermark = high_watermark
        else:
            # flush
            watermark = 1
        while len(tbd) >= watermark:
            delete_now, tbd = tbd[:batch_size], tbd[batch_size:]
            max_threads.acquire()
            threading.Thread(target=delete_many, args=[delete_now, max_threads]).start()
        to_be_deleted[0] = tbd

def async_deletes_flush():
    async_delete(None)

def get_sequences_by_accession_list_from_file(accession2seq, nt_loc_dict, nt_file):
    with open(nt_file) as ntf:
        for accession_id, accession_info in accession2seq.iteritems():
            accession_info['seq'] = get_sequence_by_accession_id_ntf(accession_id, nt_loc_dict, ntf)

def get_sequences_by_accession_list_from_s3(accession_id_list, nt_loc_dict, nt_s3_path, num_slaves, slave_id):
    threads = []
    error_flags = {}
    semaphore = threading.Semaphore(8)
    mutex = threading.RLock()
    nt_bucket, nt_key = nt_s3_path[5:].split("/", 1)
    for accession_id, accession_info in accession_id_list.iteritems():
        semaphore.acquire()
        t = threading.Thread(
            target=get_sequence_for_thread,
            args=[error_flags, accession_info, accession_id, nt_loc_dict, nt_bucket, nt_key, semaphore, mutex, async_delete, num_slaves, slave_id]
        )
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    async_deletes_flush()
    if error_flags:
        raise "Sorry there was an error"

def get_sequence_for_thread(error_flags, accession_info, accession_id, nt_loc_dict, nt_bucket, nt_key, semaphore, mutex, async_delete_func, num_slaves, slave_id, seq_count=[0]): #pylint: disable=dangerous-default-value
    try:
        seq = get_sequence_by_accession_id_s3(accession_id, nt_loc_dict, nt_bucket, nt_key, async_delete_func)
        with mutex:
            accession_info['seq'] = seq
            seq_count[0] += 1
            if slave_id == 0 and seq_count[0] % 100 == 0:
                print("%d sequences fetched by slave 0 of %d, most recently %s" % (seq_count[0], num_slaves, accession_id))
    except:
        with mutex:
            if not error_flags:
                traceback.print_exc()
            error_flags["error"] = 1
    finally:
        semaphore.release()

def get_sequence_by_accession_id_ntf(accession_id, nt_loc_dict, ntf):
    ref_seq = 'NOT FOUND'
    entry = nt_loc_dict.get(accession_id)
    if entry:
        range_start = entry[0]+entry[1]
        seq_len = entry[2]
        ntf.seek(range_start, 0)
        ref_seq = ntf.read(seq_len).replace("\n", "")
    return ref_seq

def get_sequence_by_accession_id_s3(accession_id, nt_loc_dict, nt_bucket, nt_key, delete_func=os.remove):
    ref_seq = 'NOT FOUND'
    entry = nt_loc_dict.get(accession_id)
    if entry:
        range_start = entry[0]+entry[1]
        seq_len = entry[2]
        tmp_file = '/tmp/accession-%s' % accession_id
        os.mkfifo(tmp_file)
        get_range = "aws s3api get-object --range bytes=%d-%d --bucket %s --key %s %s" % (range_start, range_start + seq_len - 1, nt_bucket, nt_key, tmp_file)
        getter_subproc = subprocess.Popen(get_range, shell=True, stdout=DEVNULL)
        with open(tmp_file, "rb") as tf:
            ref_seq = tf.read().replace("\n", "")
        exitcode = getter_subproc.wait()
        assert exitcode == 0
        delete_func(tmp_file)
    return ref_seq

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
