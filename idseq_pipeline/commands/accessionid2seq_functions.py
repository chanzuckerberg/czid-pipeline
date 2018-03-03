import subprocess
import shelve
import re
import json
import threading
import traceback
import os

from .common import * #pylint: disable=wildcard-import

REF_DISPLAY_RANGE = 100
MAX_SEQ_DISPLAY_SIZE = 6000


# Test with the following function call
# generate_alignment_viz_json('../../nt','nt.db','NT', 'taxids.gsnapl.unmapped.bowtie2.lzw.cdhitdup.priceseqfilter.unmapped.star.m8', 'taxid_annot_sorted_nt.fasta', 'align_viz')

def generate_alignment_viz_json(nt_file, nt_loc_db, db_type,
                                annotated_m8, annotated_fasta,
                                output_json_dir):
    """Generate alignment details from the reference sequence, m8 and annotated fasta """
    # Go through annotated_fasta with a db_type (NT/NR match).  Infer the family/genus/species info

    if db_type != 'NT' and db_type != 'NR':
        return

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
    print("Read to Seq dictionary size: %d" % len(read2seq))

    # Go through m8 file infer the alignment info. grab the fasta sequence, lineage info, construct a tree, get the actual sequence from nt_file
    nt_loc_dict = shelve.open(nt_loc_db)
    (ntf, nt_s3_path) = (None, None)
    if nt_file.startswith("s3://"):
        nt_s3_path = nt_file
    else:
        ntf = open(nt_file)

    result_dict = {}
    line_count = 0
    accession2seq = {}
    with open(annotated_m8, 'r') as m8f:
        for line in m8f:
            line_columns = line.rstrip().split("\t")
            accession_id = line_columns[1]
            accession2seq[accession_id] = accession2seq.get(accession_id, {'count': 0})
            accession2seq[accession_id]['count'] += 1
    print("%d unique accession ids" % len(accession2seq))
    if ntf:
        for accession_id, accession_info in accession2seq.iteritems():
            accession_info['seq'] = get_sequence_by_accession_id(accession_id, nt_loc_dict, ntf, None)
    else:
        get_sequences_by_accession_list_s3(accession2seq, nt_loc_dict, nt_s3_path)

    with open(annotated_m8, 'r') as m8f:
        for line in m8f:
            line_count += 1
            line_columns = line.rstrip().split("\t")
            read_id = extract_m8_readid(line_columns[0])
            accession_id = line_columns[1]
            metrics = line_columns[2:]
            seq_info = read2seq.get(read_id)
            if line_count % 1000 == 0:
                print("%d lines in the m8 file processed." % line_count)

            if seq_info:
                [sequence, family_id, genus_id, species_id] = seq_info
                result_dict[family_id] = result_dict.get(family_id, {})
                result_dict[family_id][genus_id] = result_dict[family_id].get(genus_id, {})
                result_dict[family_id][genus_id][species_id] = result_dict[family_id][genus_id].get(species_id, {})
                accession_dict = result_dict[family_id][genus_id][species_id].get(accession_id, {})
                ref_seq = accession2seq[accession_id]['seq']
                accession_dict['reads'] = accession_dict.get('reads', [])
                ref_start = int(metrics[-4])
                ref_end = int(metrics[-3])
                if ref_start > ref_end: # SWAP
                    (ref_start, ref_end) = (ref_end, ref_start)
                ref_start -= 1
                prev_start = (ref_start - REF_DISPLAY_RANGE) if (ref_start - REF_DISPLAY_RANGE) > 0 else 0
                post_end = ref_end + REF_DISPLAY_RANGE
                if len(ref_seq) > MAX_SEQ_DISPLAY_SIZE:
                    accession_dict['ref_seq'] = '...Reference Seq Too Long ...'
                else:
                    accession_dict['ref_seq'] = ref_seq
                ref_read = [ref_seq[prev_start:ref_start], ref_seq[ref_start:ref_end], ref_seq[ref_end:post_end]]
                accession_dict['reads'].append((read_id, sequence, metrics, ref_read))
                accession_dict['ref_seq_len'] = len(ref_seq)
                accession_dict['ref_link'] = "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=fasta" % accession_id
                result_dict[family_id][genus_id][species_id][accession_id] = accession_dict
    print("%d lines in the m8 file" % line_count)

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
    return "Read2Seq Size: %d, M8 lines %d, %d unique accession ids" % (len(read2seq), line_count, len(accession2seq))

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

def async_delete(tmp_file, mutex=threading.RLock(), to_be_deleted=[[]], batch_size=512, max_threads=threading.Semaphore(8)): #pylint: disable=dangerous-default-value
    with mutex:
        tbd = to_be_deleted[0]
        if tmp_file:
            tbd.append(tmp_file)
            watermark = batch_size
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

def get_sequences_by_accession_list_s3(accession_id_list, nt_loc_dict, nt_s3_path):
    threads = []
    error_flags = {}
    semaphore = threading.Semaphore(64)
    mutex = threading.RLock()
    for accession_id, accession_info in accession_id_list.iteritems():
        semaphore.acquire()
        t = threading.Thread(
            target=get_sequence_for_thread,
            args=[error_flags, accession_info, accession_id, nt_loc_dict, nt_s3_path, semaphore, mutex, async_delete]
        )
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    async_deletes_flush()
    if error_flags:
        raise "Sorry there was an error"

def get_sequence_for_thread(error_flags, accession_info, accession_id, nt_loc_dict, nt_s3_path, semaphore, mutex, async_delete_func, seq_count=[0]): #pylint: disable=dangerous-default-value
    try:
        seq = get_sequence_by_accession_id(accession_id, nt_loc_dict, None, nt_s3_path, async_delete_func)
        with mutex:
            accession_info['seq'] = seq
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

def get_sequence_by_accession_id(accession_id, nt_loc_dict, ntf, nt_s3_path, async_delete_func=None):
    if not ntf: # Had to open individual connections to be thread safe
        (nt_bucket, nt_key) = nt_s3_path[5:].split("/", 1)
    entry = nt_loc_dict.get(accession_id)
    if entry:
        range_start = entry[0]+entry[1]
        seq_len = entry[2]
        if ntf:
            ntf.seek(range_start, 0)
            ref_seq = ntf.read(seq_len)
        else:
            tmp_file = 'accession-%s' % accession_id
            get_range = "aws s3api get-object --range bytes=%d-%d --bucket %s --key %s %s" % (range_start, range_start + seq_len - 1, nt_bucket, nt_key, tmp_file)
            subprocess.check_output(get_range, shell=True)
            with open(tmp_file) as af:
                ref_seq = af.read()
            if async_delete_func:
                async_delete_func(tmp_file)
            else:
                os.remove(tmp_file)
            # use AWS api
        return ref_seq.replace("\n", "")
    else:
        return 'NOT FOUND'

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
    summary = generate_alignment_viz_json(db_path, local_db_loc_path, db_type,
                                          local_m8_path, local_fasta_path, local_json_path)
    summary_file_name = "%s.summary" % local_json_path
    with open(summary_file_name, 'w') as summaryf:
        summaryf.write(summary)
    # copy the data over
    execute_command("aws s3 cp --quiet %s %s --recursive" % (local_json_path, output_json_s3_path))
    execute_command("aws s3 cp --quiet %s %s/" % (summary_file_name, os.path.dirname(output_json_s3_path)))

    # Clean up
    execute_command("rm -rf %s" % dest_dir)
