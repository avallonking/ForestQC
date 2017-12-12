import numpy as np
import pandas as pd
import statistics
from operator import itemgetter

DP_GQ_START_IDX = 9
GENOTYPE_IDX = 0
GENOTYPE_DEPTH_IDX = 1
DP_IDX = 2
GQ_IDX = 3
# DP_THRESHOLD = 34
# GQ_THRESHOLD = 99

# for PED file processing
FAMILY_ID_IDX = 0
INDIVIDUAL_ID_IDX = 1
FATHER_ID_IDX = 2
MOTHER_ID_IDX = 3
SEX_IDX = 4
PHENOTYPE_IDX = 5
SAMPLE_ID_IDX = -1

def isHomozygous(genotype):
    if genotype[0] == genotype[-1]:
        return True
    else:
        return False

def statDPGQ(variant_info, dp_threshold=34, gq_threshold=99, target_idx=None, chr=None):
    # get Mean_DP, Mean_GQ, SD_DP, SD_GQ, Outlier_DP, Outlier_GQ from a line of variant information
    dp = []
    gq = []
    outlier_dp_count = 0
    outlier_gq_count = 0

    individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]
    if target_idx and chr == 'chrX':
      individual_info = itemgetter(*target_idx)(individual_info) 

    for idv in individual_info:
        idv_info = idv.split(':')
        if len(idv_info) < 4:
          continue

        if idv_info[DP_IDX] != '.':
          dp.append(float(idv_info[DP_IDX]))
          if float(idv_info[DP_IDX]) < dp_threshold:
            outlier_dp_count += 1
          # if idv_info[DP_IDX] != '0':
            # avg_dp.append(float(idv_info[DP_IDX]))

        if idv_info[GQ_IDX] != '.':
          gq.append(float(idv_info[GQ_IDX]))
          if float(idv_info[GQ_IDX]) < gq_threshold:
            outlier_gq_count += 1
          # if idv_info[GQ_IDX] != '0':
            # avg_gq.append(float(idv_info[GQ_IDX]))

    if len(gq) > 1:
      mean_gq = round(statistics.mean(gq), 5)
      sd_gq = round(statistics.stdev(gq), 5)
      outlier_gq = round((outlier_gq_count / len(gq)), 5)
    elif len(gq) == 1:
      mean_gq = gq[0]
      sd_gq = 0
      outlier_gq = outlier_gq_count
    else:
      mean_gq = 'NA'
      sd_gq = 'NA'
      outlier_gq = 'NA'

    if len(dp) > 1:
      mean_dp = round(statistics.mean(dp), 5)
      sd_dp = round(statistics.stdev(dp), 5)
      outlier_dp = round((outlier_dp_count) / len(dp), 5)
    elif len(dp) == 1:
      mean_dp = dp[0]
      sd_dp = 0
      outlier_dp = outlier_dp_count
    else:
      mean_dp = 'NA'
      sd_dp = 'NA'
      outlier_dp = 'NA'

    return mean_dp, mean_gq, sd_dp, sd_gq, outlier_dp, outlier_gq

def getGC_table(gc_file):
  gc_table = pd.read_table(gc_file, header=None)
  gc_table.columns = ['CHR', 'Start_Pos', 'GC']
  chr_gc = {}
  for i in set(gc_table.CHR):
    chr_gc[i] = gc_table[gc_table.CHR == i]
  return chr_gc

def getGC(pos, gc_table, step=1000):
  gc = gc_table.iloc[int(int(pos) / step), 2]
  return gc

def computeAB(ref, alt):
    # calculate ABHom or ABHet
    if ref + alt == 0:
        return 'NA'
    else:
        return round(ref / (ref + alt), 5)

def getSexInfo(ped_file, gender_file):
  # read a ped file and return a list of male sample ids and female sample ids
  # assume we only have either ped_file or gender_file as input
  male = []
  female = []

  try:
    with open(gender_file, 'r') as g:
      next(g)
      for line in g:
        info = line.strip().split('\t')
        sample_id = info[0]
        gender = info[1]
        if gender == 'm':
          male.append(sample_id)
        elif gender == 'f':
          female.append(sample_id)
  except TypeError:
    pass

  try:
    with open(ped_file, 'r') as p:
      for line in p:
        info = line.strip().split('\t')
        sample_id = info[SAMPLE_ID_IDX]
        if sample_id == 'NA':
          continue
        sex = info[SEX_IDX]
        if sex == 'm':
          male.append(sample_id)
        elif sex == 'f':
          female.append(sample_id)
  except TypeError:
    pass

  return male, female

def getTargetIdx(male, female, sample_list):
  # male, female = getSexInfo(ped_file)
  male_idx = list(map(sample_list.index, male))
  female_idx = list(map(sample_list.index, female))
  return male_idx, female_idx

def getAB(variant_info, sample_level_AB=False, chr=None, target_idx=None):
    # get ABHet and ABHom from a line of variant information
    sampleAB = {}
    variantAB = {'abhet': [0, 0], 'abhom': [0, 0]}
    sample_idx = 1

    individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]
    if target_idx and chr == 'chrX':
      individual_info = itemgetter(*target_idx)(individual_info) 

    for idv in individual_info:
        genotype = idv.split(':')[GENOTYPE_IDX]

        if '.' in genotype:
            continue

        allele_depth = idv.split(':')[GENOTYPE_DEPTH_IDX].split(',')

        if isHomozygous(genotype):
            sample_info = 'ABHom' + str(sample_idx)
            sample_idx += 1

            variantAB['abhom'][0] += float(allele_depth[0])
            variantAB['abhom'][1] += float(allele_depth[-1])
        else:
            sample_info = 'ABHet' + str(sample_idx)
            sample_idx += 1

            variantAB['abhet'][0] += float(allele_depth[0])
            variantAB['abhet'][1] += float(allele_depth[-1])

        if sample_level_AB:
            sampleAB[sample_info] = computeAB(float(allele_depth[0]), float(allele_depth[1]))

    ABHet = computeAB(variantAB['abhet'][0], variantAB['abhet'][1])
    ABHom = computeAB(variantAB['abhom'][0], variantAB['abhom'][1])

    if sample_level_AB:
        return ABHet, ABHom, sampleAB
    else:
        return ABHet, ABHom


def getMAF(variant_info, target_idx=None, chr=None):
    allele_info = {}
    freq = []

    individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]
    if target_idx and chr == 'chrX':
      individual_info = itemgetter(*target_idx)(individual_info) 

    for idv in individual_info:
        genotype = idv.split(':')[GENOTYPE_IDX]
        for i in [0, -1]:
            if genotype[i] == '.':
                continue
            if genotype[i] not in allele_info:
                allele_info[genotype[i]] = 1
            else:
                allele_info[genotype[i]] += 1
    for info in allele_info.values():
        freq.append(info / sum(allele_info.values()))
        
    try:
      maf = sorted(freq)[-2] if len(freq) >= 2 else freq[-1]
      maf = round(maf, 5)
    except IndexError:
      maf = 'NA'
    return maf


def getMissing(variant_info, chr=None, male_idx=None, target_idx=None):
    # get missing rate from a line of variant
    missing_allele = 0

    individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]
    if target_idx and chr == 'chrX':
      individual_info = itemgetter(*target_idx)(individual_info) 

    total_allele = len(individual_info) * 2
    for idv in individual_info:
        genotype = idv.split(':')[GENOTYPE_IDX]

        allele_first = genotype[0]
        allele_second = genotype[-1]
        if allele_first == '.':
            missing_allele += 1
        if allele_second == '.':
            missing_allele += 1

    if chr == 'chrX' and male_idx:
      male_info = itemgetter(*male_idx)(individual_info)
      for idv in male_info:
        genotype = idv.split(':')[GENOTYPE_IDX]
        if '.' not in genotype and not isHomozygous(genotype):
          missing_allele += 1

    missing_rate = round(missing_allele / total_allele, 5)
    return missing_rate


def getHWE_Direct(hwe_file):
  # get HWE p-value directly from the file provided by users
  # input: a file have 2 tab-separated columns: SNP_ID and HWE_p-value
  # output: a dictionary contianing the information
  hwe_dict = {}
  try:
    with open(hwe_file, 'r') as h:
      for line in h:
        info = line.strip().split('\t')
        snp = info[0]
        hwe = info[1]
        try:
          hwe_dict[snp] = float(hwe)
        except ValueError:
          continue
  except TypeError:
    pass
  return hwe_dict

def getControlSamples(ped_file, sample_list):
  # get sample information (control, case or missing) from pedigree file
  # we only need the control samples
  # output: the index of control samples in sample list
  control_samples = []
  control_samples_idx = []
  try:
    with open(ped_file, 'r') as p:
      for line in p:
        info = line.strip().split('\t')
        sample_id = info[SAMPLE_ID_IDX]
        phenotype_id = info[PHENOTYPE_IDX]
        if phenotype_id == '1' and sample_id != 'NA':
          control_samples.append(sample_id)
    
    control_samples_idx = list(map(sample_list.index, control_samples))
  except TypeError:
    pass
  
  return control_samples_idx

def getHWE(variant_info, control_samples_idx=None, chr=None, target_idx=None):
    # get Hardy-Weinberg Equilibrium P-value from a line of variant
    hom0 = 0  # homozygous
    hom1 = 0
    het = 0  # heterozygous
    total_non_missing_haploid = 0
    het_probs = []
    hwe = 0.0

    individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]

    if control_samples_idx:
      female_control = list(set(control_samples_idx).intersection(target_idx))

    if chr == 'chrX':
      if control_samples_idx:
        individual_info = itemgetter(*female_control)(individual_info)
      elif target_idx:
        individual_info = itemgetter(*target_idx)(individual_info)
    elif control_samples_idx:
      individual_info = itemgetter(*control_samples_idx)(individual_info)

    for idv in individual_info:
        genotype = idv.split(':')[GENOTYPE_IDX]
        if '.' not in genotype:
            total_non_missing_haploid += 1
            if genotype == '0/0' or genotype == '0|0':
                hom0 += 1
            elif genotype == '1/1' or genotype == '1|1':
                hom1 += 1
            else:
                het += 1
    hom_rare = hom0 if hom0 < hom1 else hom0
    rare_copies = 2 * hom_rare + het
    mid = int(rare_copies * (2 * total_non_missing_haploid - rare_copies) / (2 * total_non_missing_haploid))
    for i in range(rare_copies + 1):
        het_probs.append(0.0)

    if ((rare_copies & 1) ^ (mid & 1)):
        mid += 1

    het_probs[mid] = 1.0
    prob_sum = het_probs[mid]

    curr_het = mid
    curr_hom_rare = (rare_copies - mid) / 2
    curr_hom_common = total_non_missing_haploid - curr_het - curr_hom_rare
    for curr_het in range(mid, 1, -2):
        het_probs[curr_het - 2] = het_probs[curr_het] * curr_het * (curr_het - 1.0) / (
        4.0 * (curr_hom_rare + 1.0) * (curr_hom_common + 1.0))
        prob_sum += het_probs[curr_het - 2]
        curr_hom_common += 1
        curr_hom_rare += 1

    curr_het = mid
    curr_hom_rare = (rare_copies - mid) / 2
    curr_hom_common = total_non_missing_haploid - curr_het - curr_hom_rare
    for curr_het in range(mid, rare_copies, 2):
        if curr_het > rare_copies - 2:
            continue
        het_probs[curr_het + 2] = het_probs[curr_het] * 4.0 * curr_hom_rare * curr_hom_common / (
        (curr_het + 2.0) * (curr_het + 1.0))
        prob_sum += het_probs[curr_het + 2]
        curr_hom_rare -= 1
        curr_hom_common -= 1

    for t in range(rare_copies + 1):
        het_probs[t] /= prob_sum

    for r in range(rare_copies + 1):
        if het_probs[r] > het_probs[het]:
            continue
        hwe += het_probs[r]
    hwe = 1.0 if hwe > 1 else round(hwe, 5)

    return hwe


def getFamilyRelation(ped_file, sample_list):
    # get family relation from ped file, return a dict that contains each individual and his relationships with other samples
    # intermediate format: { FamilyID: { individualID: [ FatherID, MotherID, seqID ] } }
    # final output format: { individual_sample_id: [ father_sample_id, mother_sample_id ] }
    pedigree = {}
    relationship = {}
    try:
      with open(ped_file, 'r') as ped:
          for line in ped:
              if line.startswith('Family'):
                  continue
              id_set = line.strip().split('\t')
              family_id = id_set[FAMILY_ID_IDX]
              individual_id = id_set[INDIVIDUAL_ID_IDX]
              father_id = id_set[FATHER_ID_IDX]
              mother_id = id_set[MOTHER_ID_IDX]
              sample_id = id_set[SAMPLE_ID_IDX] if id_set[SAMPLE_ID_IDX] in sample_list else 'NA'
              if family_id not in pedigree:
                  pedigree[family_id] = {}
              pedigree[family_id][individual_id] = [father_id, mother_id, sample_id]
      for member in pedigree.values():
          for individual_relation in member.values():
              individual_sample_id = individual_relation[-1]
              father_sample_id = 'NA' if individual_relation[0] not in member else member[individual_relation[0]][-1]
              mother_sample_id = 'NA' if individual_relation[1] not in member else member[individual_relation[1]][-1]
              if 'NA' in [father_sample_id, mother_sample_id] or individual_sample_id == 'NA':
                  continue
              relationship[individual_sample_id] = [father_sample_id, mother_sample_id]
    except TypeError:
      pass
    return relationship

def getMendel(variant_info, sample_list, relationship, chr=None, male_list=None):
    # get mendel error rate of a variant from one line of variant information
    # sample_list is a list of sample names
    if not relationship:
      return 'NA'
    if not male_list:
      male_list = []

    genotype_set = []
    mendel_error = 0
    individual_info = variant_info.split('\t')[DP_GQ_START_IDX:]
    sample_size = len(individual_info)

    for idv in individual_info:
        genotype = idv.split(':')[GENOTYPE_IDX]
        genotype_set.append(genotype)
    sample_genotype = dict(zip(sample_list, genotype_set))

    for child, parents in relationship.items():
        child_genotype = sample_genotype[child]
        father_genotype = sample_genotype[parents[0]]
        mother_genotype = sample_genotype[parents[1]]
        if chr == 'chrX':
          if child in male_list:
            # if '.' in child_genotype +  mother_genotype:
              # continue
            # if isHomozygous(child_genotype) and child_genotype[0] not in mother_genotype:
              # mendel_error += 1
            continue
          else:
            if '.' in child_genotype + father_genotype + mother_genotype:
                continue
            normal_inherited_genotype = [father_genotype[:2] + mother_genotype[0],
                                         father_genotype[:2] + mother_genotype[-1],
                                         father_genotype[:0:-1] + mother_genotype[0],
                                         father_genotype[:0:-1] + mother_genotype[-1]]
            if child_genotype not in normal_inherited_genotype and child_genotype[::-1] not in normal_inherited_genotype:
              mendel_error += 1
        else:
          if '.' in child_genotype + father_genotype + mother_genotype:
              continue
          normal_inherited_genotype = [father_genotype[:2] + mother_genotype[0],
                                       father_genotype[:2] + mother_genotype[-1],
                                       father_genotype[:0:-1] + mother_genotype[0],
                                       father_genotype[:0:-1] + mother_genotype[-1]]
          if child_genotype not in normal_inherited_genotype and child_genotype[::-1] not in normal_inherited_genotype:
            mendel_error += 1
    mendel_error_rate = round(mendel_error / sample_size, 5)
    return mendel_error_rate

def getDiscordantGenotype(variant_info, discordant_genotype_dict):
    rsid = variant_info.split('\t')[0] + ':' + variant_info.split('\t')[1]
    snp_id = variant_info.split('\t')[2]
    try:
        return discordant_genotype_dict[rsid]
    except KeyError:
        pass

    try:
        return discordant_genotype_dict[snp_id]
    except KeyError:
        return 'NA'

def get_additional_features(features_df, rsid):
    # this input is a pandas dataframe
    try:
        return list(features_df[features_df['RSID'] == rsid].iloc[0, 1:])
    except IndexError:
        return ['NA'] * (features_df.shape[1] - 1)
    except TypeError:
        pass
