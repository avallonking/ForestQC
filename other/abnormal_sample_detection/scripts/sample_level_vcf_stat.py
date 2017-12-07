import os
import glob
import gzip
import numpy as np

def isNotMissing(genotype):
  if '.' in genotype:
    return False
  else:
    return True

def isHeterozygous(genotype):
  if genotype[0] != genotype[-1]:
    return True
  else:
    return False

def selectedVariants(chr, pos, list_):
  # input: chromosome, position and list_ that contains needed variants (RSID)
  #        chr and pos are strings
  #        list_ is a list
  # return: True if this variant is in the list_, else False
  if 'chr' not in chr.lower():
    chr = 'chr' + chr
  rsid = chr + ':' + pos

  try:
    if rsid in list_:
      return True
    else:
      return False
  except TypeError:
    return True

def getSampleIdx(info_line, sample):
  sample_idx = {}
  sample_list = info_line.strip().split('\t')[9:]
  for individual in sample:
    try:
      sample_idx[individual] = sample_list.index(individual)
    except ValueError:
      continue
  return sample_idx

def sampleLevelAB(vcf_list, subset = None):
  # only ABHet
  sample_allele = {}
  # ab = {}

  for vcf in vcf_list:
    with gzip.open(vcf, 'rt') as v:
      for line in v:
        if line.startswith('#CHROM'):
          sample = line.strip().split('\t')[9:]
          for i in sample:
            sample_allele[i] = {'ALT': 0, 'REF': 0}

        elif not line.startswith('#'):
          entry = line.strip().split('\t')
          info = entry[9:]
          for individual, sample_info in zip(sample, info):
            genotype = sample_info.split(':')[0]
            depth = sample_info.split(':')[1]
            
            if isNotMissing(genotype) and isHeterozygous(genotype):
              sample_allele[individual]['REF'] += float(depth.split(',')[0])
              sample_allele[individual]['ALT'] += float(depth.split(',')[1])

  return sample_allele

def sample_ab_calculation(sample_allele):
  ab = {}
  for k, v in sample_allele.items():
    try:
      ab[k] = round(v['REF'] / (v['ALT'] + v['REF']), 5)
    except ZeroDivisionError:
      continue

  return ab

def sampleMissing(vcf_list, sample, list_=None):
  missing = {}

  for vcf in vcf_list:
    with gzip.open(vcf, 'rt') as v:
      for line in v:
        if line.startswith('#CHROM'):
          sample_idx = getSampleIdx(line, sample)

        elif line.startswith('chr'):
          entry = line.strip().split('\t')
          info = entry[9:]
          if not selectedVariants(entry[0], entry[1], list_):
            continue
          for individual2, individual2_idx in sample_idx.items():
            sample_info = info[individual2_idx]
            genotype = sample_info.split(':')[0]
            missing_allele = genotype.count('.')
            try:
              missing[individual2]['missing_allele'] += missing_allele
              missing[individual2]['total_allele'] += 1
            except KeyError:
              missing[individual2] = {'missing_allele' : missing_allele, 'total_allele': 0}

  return missing

def sampleDP(vcf_list, sample):
  # get mean and SD of depth of each individual
  dp = {}
  for vcf in vcf_list:
    with gzip.open(vcf, 'rt') as v:
      for line in v:
        if line.startswith('#CHROM'):
          sample_idx = getSampleIdx(line, sample)

        elif line.startswith('chr'):
          entry = line.strip().split('\t')
          info = entry[9:]
          for individual, individual_idx in sample_idx.items():
            sample_info = info[individual_idx]
            try:
              sample_dp = float(sample_info.split(':')[2])
            except ValueError:
              sample_dp = np.nan
            except IndexError:
              allele_dp = sample_info.split(':')[1]
              sample_dp = float(allele_dp.split(',')[0]) + float(allele_dp.split(',')[1])

            try:
              dp[individual].append(sample_dp)
            except KeyError:
              dp[individual] = [sample_dp]
  return dp

def scanVCF(dir, sample, subset=None):
  # sample should be a list, dir should be a string
  file_list = []
  chr_list = os.listdir(dir)
  for chr in chr_list:
    if chr == 'chrX':
      continue
    temp = glob.glob(os.path.join(dir, chr, 'orig', '*.vcf.gz'))
    temp2 = [elem for elem in temp if 'site_info' not in elem]
    file_list.extend(temp2)
  return sampleLevelAB(vcf_list=file_list, sample=sample, subset=subset)
