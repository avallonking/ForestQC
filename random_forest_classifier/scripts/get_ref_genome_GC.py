# get GC content with 1000bp sliding windows.
# usage:
#   python3 get_ref_genome_GC.py ref_genome.fasta output_filename
# input: 
#   ref: reference genome in fasta format
#   out: output file name
# output:
#   a dataframe, which has CHR, Start_Position and GC columns

import sys
import numpy as np

def computeGC(seq):
  try:
    gc = (seq.count('G') + seq.count('C')) / (seq.count('G') + seq.count('C') + seq.count('A') + seq.count('T'))
    return gc
  except ZeroDivisionError:
    return np.nan

def main():
  # get the reference genome
  ref = sys.argv[1]
  out = sys.argv[2]
  r = open(ref,'r')
  o = open(out, 'w')

  seq_list = []
  start_pos = 1
  for line in r:
    if line.startswith('>'):
      if seq_list:
        seq = ''.join(seq_list)
        gc = computeGC(seq)
        seq_list = []
        o.write('\t'.join([chr,str(start_pos),str(gc)]) + '\n')
      chr = line.strip()[1:]
      start_pos = 1
    else:
      seq_list.append(line.strip().upper())
      if len(seq_list) == 20:
        seq = ''.join(seq_list)
        gc = computeGC(seq)
        seq_list = []
        o.write('\t'.join([chr,str(start_pos),str(gc)]) + '\n')
        start_pos += 1000
  r.close()
  o.close()

if __name__ == '__main__':
  main()
