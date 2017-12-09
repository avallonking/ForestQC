# This is tool is deprecated

import sys
import gzip

def indelRatio(vcf):
  insertion = 0
  deletion = 0

  if type(vcf) is list: 
    for file in vcf:
      with gzip.open(file, 'rt') as v:
        for line in v:
          if not line.startswith('#'):
            info = line.strip().split('\t')
            ref = info[3]
            alt = info[4]
            if len(ref) <= len(alt):
              # insertion
              insertion += 1
            else:
              deletion += 1

  else:
    with open(vcf, 'r') as v:
      for line in v:
        if not line.startswith('#'):
          info = line.strip().split('\t')
          ref = info[3]
          alt = info[4]
          if len(ref) <= len(alt):
            # insertion
            insertion += 1
          else:
            deletion += 1
  return insertion / deletion

def main():
  vcf = sys.argv[1]
  print('Indel Ratio: ' + str(indelRatio(vcf)))

if __name__ == '__main__':
  main()
