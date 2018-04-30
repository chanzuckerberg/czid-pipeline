from common import *
import shelve

def remove_short(input_fasta, length_threshold, size_db_file, size_count_file, output_fasta, removed_fasta):
    '''Remove fasta records with sequence shorter than threshold; report length distribution'''
    size_count = shelve.open(size_db_file)
    size_count.clear()
    with open(input_fasta, 'rb') as input_f:
        with open(output_fasta, 'wb') as output_f:
            with open(removed_fasta, 'wb') as removed_f:
                line_number = 0
                def write_output(seq_name, seq, seq_len):
                    size_count[str(seq_len)] = size_count.get(str(seq_len), 0) + 1
                    if seq_len >= length_threshold:
                        output_f.write(seq_name + seq)
                    else:
                        removed_f.write(seq_name + seq)
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
                        seq_len += len(line.replace("\n",""))
                write_output(seq_name, seq, seq_len)
    with open(size_count_file, 'wb') as count_f:
        for seqlen in sorted(size_count.keys(), key=lambda x: int(x)):
            count_f.write("sequence length %s: %d records" % (seqlen, size_count[seqlen]))
    size_count.close()

def main():
    dest_dir = os.path.join(DEST_DIR, 'curate_nt')
    execute_command("mkdir -p %s" % dest_dir)

    nt_s3 = "s3://idseq-database/alignment_data/2018-04-01-utc-1522569777-unixtime__2018-04-04-utc-1522862260-unixtime/nt"
    nt_local = "%s/nt" % dest_dir
    if not os.path.isfile(nt_local):
        execute_command("s3mi cp %s %s" % (nt_s3, nt_local))

    length_threshold = 100
    size_db_file = "%s/nt_seqlen.db" % dest_dir
    size_count_file = "%s/nt_length_counts.txt" % dest_dir
    output_fasta = "%s/nt_curated.fasta" % dest_dir
    removed_fasta = "%s/nt_removed.fasta" % dest_dir
    
    remove_short(nt_local, length_threshold, size_db_file, size_count_file, output_fasta, removed_fasta)

if __name__ == "__main__":
    main()
