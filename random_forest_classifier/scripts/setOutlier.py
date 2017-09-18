import sys
import gzip
import numpy as np
import random
from operator import itemgetter

DP_GQ_START_IDX = 9
DP_IDX = 2
GQ_IDX = 3
DP_THRESHOLD = 0
GQ_THRESHOLD = 0

def main():
  file_list = sys.argv[1:]
  gq = []
  dp = []
  # sample_idx = None
  for file in file_list:
    f = gzip.open(file, 'rt') if 'gz' in file else open(file, 'r')
    for line in f:
      if not line.startswith('#'):
        variant_info = line
        individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]
        # if not sample_idx:
          # sample_idx = random.sample(range(len(individual_info)), int(len(individual_info) / 10))
        # individual_info = itemgetter(*sample_idx)(individual_info)
        for idv in individual_info:
          idv_info = idv.split(':')
          if len(idv_info) < 4:
            continue
          if idv_info[DP_IDX] != '.':
            dp.append(float(idv_info[DP_IDX]))
          if idv_info[GQ_IDX] != '.':
            gq.append(float(idv_info[GQ_IDX]))

    f.close()

  DP_THRESHOLD = np.percentile(dp, 25, interpolation = 'lower')
  GQ_THRESHOLD = np.percentile(gq, 25, interpolation = 'lower')
  print('Outlier_DP: ' + str(DP_THRESHOLD))
  print('Outlier_GQ: ' + str(GQ_THRESHOLD))

if __name__ == '__main__':
  main()
