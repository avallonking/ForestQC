import sys
import gzip
import pandas as pd
import numpy as np

def getDP(sample_info):
  try:
    dp = float(sample_info.split(':')[2])
  except:
    dp = np.nan

  return dp

def main():
  region = 0
  temp = {}
  result = {}
  region_size = {}
  input = sys.argv[1]
  output = sys.argv[2]
  step = int(sys.argv[3])

  # with open('/u/scratch2/k/k8688933/dp_stat_vcf_only/dp_csv_split2/header', 'r') as h:
    # header = h.readline().strip().split(',')


  with gzip.open(input, 'rt') as f:
    for line in f:

      if line.startswith('#CHROM'):
        header = line.strip().split('\t')[9:]
        for sample in header:
          result[sample] = {}

      elif not line.startswith('#'):
        info_ = line.strip().split('\t')
        sample_info_list = info_[9:]
        dp = list(map(getDP, sample_info_list))
        pos = float(info_[1])
        current_region = int(pos / step)

        if current_region > region:
          for k, v in temp.items():
            result[k][region] = np.nanmean(v)
            region_size[region] = len(v)

          temp = {}
          for sample, depth in zip(header, dp):

            # if depth >= 100:
              # continue

            temp[sample] = [depth]
          region = current_region

        elif current_region == region:
          for sample, depth in zip(header, dp):

            # if depth >= 100:
              # continue

            try:
              temp[sample].append(depth)
            except KeyError:
              temp[sample] = [depth]

  res = pd.DataFrame(result)
  region_info = pd.Series(region_size)
  res['region_size'] = region_info
  res.to_csv(output)

if __name__ == '__main__':
  main()
