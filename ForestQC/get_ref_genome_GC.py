# get GC content with default 1000bp sliding windows.
# usage:
#   python3 get_ref_genome_GC.py ref_genome.fasta output_filename
# input: 
#   ref: reference genome in fasta format
#   out: output file name
# output:
#   a dataframe, which has CHR, Start_Position and GC columns

def computeGC(seq):
    try:
        gc = (seq.count('G') + seq.count('C')) / (seq.count('G') + seq.count('C') + seq.count('A') + seq.count('T'))
        return round(gc, 5)
    except ZeroDivisionError:
        return 'NA'


def execute_compute_gc(ref, out, window_size):
    # get the reference genome
    print('Computing')
    r = open(ref, 'r')
    o = open(out, 'w')

    seq = ''
    start_pos = 1
    for line in r:
        if line.startswith('>'):
            if len(seq) != 0:
                gc = computeGC(seq)
                seq = ''
                o.write('\t'.join([chr, str(start_pos), str(gc)]) + '\n')
            chr = line.split()[0][1:]
            if 'chr' not in chr:
                chr = 'chr' + chr
            start_pos = 1
        else:
            seq += line.strip().upper()
            if len(seq) >= window_size:
                gc = computeGC(seq)
                seq = ''
                o.write('\t'.join([chr, str(start_pos), str(gc)]) + '\n')
                start_pos += window_size
    r.close()
    o.close()
    print('Done')
