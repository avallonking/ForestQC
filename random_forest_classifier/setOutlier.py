import sys
import gzip
import numpy as np

DP_GQ_START_IDX = 9
DP_IDX = 2
GQ_IDX = 3
DP_THRESHOLD = 0
GQ_THRESHOLD = 0

def main():
  file_list = sys.argv[1:]
  gq = []
  dp = []
  for file in file_list:
    f = gzip.open(file, 'rt') if 'gz' in file else open(file, 'r')
    for line in f:
      if not line.startswith('#'):
        variant_info = line
        individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]
        for idv in individual_info:
          idv_info = idv.split(':')
          if len(idv_info) < 4:
            continue
          if idv_info[DP_IDX] != '.' and idv_info[DP_IDX] != '0':
            dp.append(float(idv_info[DP_IDX]))
          if idv_info[GQ_IDX] != '.' and idv_info[GQ_IDX] != '0':
            gq.append(float(idv_info[GQ_IDX]))
  DP_THRESHOLD = np.percentile(dp, 25, interpolation = 'lower')
  GQ_THRESHOLD = np.percentile(gq, 25, interpolation = 'lower')
  print('Outlier_DP: ' + str(DP_THRESHOLD))
  print('Outlier_GQ: ' + str(GQ_THRESHOLD))

if __name__ == '__main__':
  main()
