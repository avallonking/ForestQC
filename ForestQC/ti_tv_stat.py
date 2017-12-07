# only use with a whole vcf file
# input file can be vcf or vcf.gz files

import sys
import gzip

def getTiTv(vcf_file):
  ti_count = 0.0
  tv_count = 0.0
  ti = ['AG', 'GA', 'CT', 'TC']
  tv = ['AC', 'CA', 'GT', 'TG', 'AT', 'TA', 'CG', 'GC']
  f = open(vcf_file, 'rt') if 'gz' not in vcf_file else gzip.open(vcf_file, 'r')
  for line in f:
    if line.startswith('#'):
      continue
    ref = line.split('\t')[3]
    alt = line.split('\t')[4]

    if len(alt) == 1:
      conversion = ref + alt
      if conversion in ti:
        ti_count += 1
      elif conversion in tv:
        tv_count += 1
    else:
      for i in alt.split(','):
        conversion = ref + i
        if conversion in ti:
          ti_count += 1
        elif conversion in tv:
          tv_count += 1
  f.close()

  return round(ti_count / tv_count, 5)

def main():
  input_file = sys.argv[1]
  ti_tv = getTiTv(input_file)
  print('Ti / Tv: ' + str(ti_tv))

if __name__ == '__main__':
  main()
