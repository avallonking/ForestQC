from sample_level_vcf_stat import *
import os
import sys
import pandas as pd

# sample = os.listdir('/u/home/k/k8688933/Jaehoon/data')
# good_variants_rsid_file = '/u/scratch2/k/k8688933/stat_output/vqsr_qc4/good.all.clfB.rsid'
vcf_file = sys.argv[1]
outfile = sys.argv[2]
# list_ = []

# with open(good_variants_rsid_file, 'r') as f:
  # for line in f:
    # if not line.startswith('RSID'):
      # list_.append(line.strip())

sample_ab = sampleLevelAB([vcf_file])
ab = pd.DataFrame(sample_ab)
ab.to_csv(outfile)
