from common import *
import shelve

def size_select(nt_local, short_threshold, long_threshold,
                size_db_file, size_count_file,
                short_records_fasta, medium_records_fasta, long_records_fasta):

def size_select(input_fasta, length_threshold, size_db_file, size_count_file, output_fasta, removed_fasta):
    '''Get sequence length distribution; split fasta based on length'''
    size_count = shelve.open(size_db_file)
    size_count.clear()
    with open(input_fasta, 'rb') as input_f:
        with open(short_records_fasta, 'wb') as short_f:
            with open(medium_records_fasta, 'wb') as medium_f:
                with open(long_records_fasta, 'wb') as long_f:
                    line_number = 0
                    def write_output(seq_name, seq, seq_len):
                        size_count[str(seq_len)] = size_count.get(str(seq_len), 0) + 1
                        record = seq_name + seq
                        if seq_len < short_threshold:
                            short_f.write(record)
                        elif seq_len < long_threshold:
                            medium_f.write(record)
                        else:
                            long_f.write(record)
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
            count_f.write("sequence length %s: %d records\n" % (seqlen, size_count[seqlen]))
    size_count.close()

def main():
    dest_dir = os.path.join(DEST_DIR, 'curate_nt')
    execute_command("mkdir -p %s" % dest_dir)

    nt_s3 = "s3://idseq-database/alignment_data/2018-04-01-utc-1522569777-unixtime__2018-04-04-utc-1522862260-unixtime/nt"
    nt_local = "%s/nt" % dest_dir
    if not os.path.isfile(nt_local):
        execute_command("s3mi cp %s %s" % (nt_s3, nt_local))

    short_threshold = 100
    long_threshold = 10**5
    size_db_file = "%s/nt_seqlen.db" % dest_dir
    size_count_file = "%s/nt_length_counts.txt" % dest_dir
    short_records_fasta = "%s/nt_short.fasta" % dest_dir
    medium_records_fasta = "%s/nt_medium.fasta" % dest_dir
    long_records_fasta = "%s/nt_long.fasta" % dest_dir

    for output in [size_db_file, size_count_file, output_fasta, removed_fasta]:
        execute_command("rm -rf %s" % output)
    size_select(nt_local, short_threshold, long_threshold,
                size_db_file, size_count_file,
                short_records_fasta, medium_records_fasta, long_records_fasta)

if __name__ == "__main__":
    main()
