# get GC content with default 1000bp sliding windows.
# usage:
#   python3 get_ref_genome_GC.py ref_genome.fasta output_filename
# input: 
#   ref: reference genome in fasta format
#   out: output file name
# output:
#   a dataframe, which has CHR, Start_Position and GC columns

import numpy as np


def computeGC(seq):
  try:
    gc = (seq.count('G') + seq.count('C')) / (seq.count('G') + seq.count('C') + seq.count('A') + seq.count('T'))
    return gc
  except ZeroDivisionError:
    return np.nan


def execute_compute_gc(ref, out, window_size):
  # get the reference genome
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
      chr = line.strip()[1:]
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
