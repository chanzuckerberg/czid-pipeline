from common import *
import shelve

def remove_short(input_fasta, length_threshold, size_db_file, output_fasta):
    '''Remove fasta records with sequence shorter than threshold; report length distribution'''
    size_count = shelve.open(size_db_file)
    size_count.clear()
    removed = set()
    with open(input_fasta, 'rb') as input_f:
        with open(output_fasta, 'wb') as output_f:
            line_number = 0
            def write_output(seq_name, seq, seq_len):
                size_count[str(seq_len)] = size_count.get(str(seq_len), 0) + 1
                if seq_len >= length_threshold:
                    output_f.write(seq_name + seq)
                else:
                    removed.add(seq_name.rstrip("\n"))
            for line in input_f:
                line_number += 1
                if line_number % 100000 == 0:
                    print "%3.1f M lines. " %  (line_number/1000000.0)
                if line_number == 1:
                    seq_name = line
                    seq = ''
                    seq_len = 0
                elif line[0] == '>':
                    write_output(seq_name, seq, seq_len) 
                    seq_name = line
                    seq = ''
                    seq_len = 0
                else:
                    seq += line 
                    seq_len += len(seq.replace("\n",""))
            write_output(seq_name, seq, seq_len)
    for seqlen in sorted(size_count.keys(), key=lambda x: int(x)):
        print "sequence length %s: %d records" % (seqlen, size_count[seqlen])
    print "records removed: "
    print removed
    size_count.close()
    return removed

def main():
    dest_dir = os.path.join(DEST_DIR, 'curate_nt')
    execute_command("mkdir -p %s" % dest_dir)

    nt_s3 = "s3://yunfang-workdir/curate_accessionid2seq/nt.sample"
    nt_local = "%s/nt_sample.fasta" % dest_dir
    execute_command("s3mi cp %s %s" % (nt_s3, nt_local))

    length_threshold = 100
    size_db_file = "%s/nt_sample_seqlen.db" % dest_dir
    output_fasta = "%s/nt_sample_curated.fasta" % dest_dir
    
    remove_short(nt_local, length_threshold, size_db_file, output_fasta)

if __name__ == "__main__":
    main()
