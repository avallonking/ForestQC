import gzip
import argparse
from vcf_stat import *

def getDiscordInfo(discord_geno_file):
    # discordant genotype file processing
    discord_geno_dict = {}
    try:
      d = open(discord_geno_file, 'r')
      for entry in d:
          snp_id = entry.split('\t')[0]
          discord_geno_num = entry.strip().split('\t')[1]
          discord_geno_dict[snp_id] = discord_geno_num
      d.close()
    except TypeError:
      pass
    return discord_geno_dict

def vcfProcessing(vcf_file, stat_file, gc_file, ped_file, discord_geno_dict, hwe_file, gender_file):
    f = gzip.open(vcf_file, 'rt') if 'gz' == vcf_file.split('.')[-1] else open(vcf_file, 'r')
    o = open(stat_file, 'w')
    for line in f:
        if line.startswith('#CHROM'):
          sample_list = line.strip().split('\t')[DP_GQ_START_IDX:]
          break
    
    relationship = getFamilyRelation(ped_file, sample_list)
    control_samples_idx = getControlSamples(ped_file, sample_list)
    male_list, female_list = getSexInfo(ped_file, gender_file)
    male_idx, female_idx = getTargetIdx(male_list, female_list, sample_list)
    hwe_info = getHWE_Direct(hwe_file)
    gc_table_by_chr = getGC_table(gc_file)

    for line2 in f:
        if not line2.startswith('#'):
          site_info = line2.split('\t')
          chr = site_info[0] if 'chr' in site_info[0].lower() else 'chr' + site_info[0]
          pos = site_info[1]
          ref = site_info[3]
          alt = site_info[4]
          rsid = chr + ':' + pos
          maf = getMAF(line2, target_idx=female_idx, chr=chr)

          mean_dp, mean_gq, sd_dp, sd_gq, outlier_dp, outlier_gq = statDPGQ(line2, dp_threshold=dp_threshold, gq_threshold=gq_threshold, chr=chr, target_idx=female_idx)
          discordant_geno = getDiscordantGenotype(line2, discord_geno_dict)
          mendel_error = getMendel(line2, sample_list=sample_list, relationship=relationship, chr=chr, male_list=male_list)
          missing_rate = getMissing(line2, chr=chr, target_idx=female_idx)
          gc = getGC(pos, gc_table_by_chr[chr])
          if not hwe_info:
            hwe = getHWE(line2, control_samples_idx, chr=chr, target_idx=female_idx)
          else:
            try:
              hwe = hwe_info[rsid]
            except KeyError:
              hwe = 'NA'
          abhet, abhom = getAB(line2, chr=chr, target_idx=female_idx)
          o.write('\t'.join([rsid, chr, pos, ref, alt, str(maf), str(mean_dp), str(mean_gq), str(sd_dp), str(sd_gq), str(outlier_dp), str(outlier_gq), str(discordant_geno), str(mendel_error), str(missing_rate), str(hwe), str(abhet), str(abhom), str(gc)]) + '\n')
    o.close()
    f.close()

# def main():
#     # target_file: the vcf file for processing
#     # stat_file: the output filename
#     # ped_file: the pedigree file
#     # discord_geno_file: the file that has the information about discordant genotype
#     # DP_THRESHOLD: depth threshold for the calculation of Outlier_DP
#     # GQ_THRESHOLD: genotype quality threshold for the calculation of Outlier_GQ
#     parser = argparse.ArgumentParser(description='Calculate statistics for one vcf file')
#     parser.add_argument('-i', '--input', dest='target_file', help='Input vcf or vcf.gz file')
#     parser.add_argument('-o', '--output', dest='stat_file', help='Output file name')
#     parser.add_argument('-c', '--gc', dest='gc_file', help='GC content of reference genome')
#     parser.add_argument('-p', '--ped', dest='ped_file', default=None, help='Pedigree file [optional]')
#     parser.add_argument('-d', '--discord_geno_file', dest='discord_geno_file', default=None, help='Discordant genotype file [optional]')
#     parser.add_argument('-w', '--hwe', dest='hwe_file', default=None, help='HWE p-value file [optional]')
#     parser.add_argument('-g', '--gender', dest='gender_info', default=None, help='Infomation about gender of each individual [optional]')
#     parser.add_argument('--dp', dest='dp', default=34, type=float, help='Depth threshold, default = 34')
#     parser.add_argument('--gq', dest='gq', default=99, type=float, help='Genotype quality threshold, default = 99')
#     args = parser.parse_args()
#
#     target_file = args.target_file
#     stat_file = args.stat_file
#     gc_file = args.gc_file
#     ped_file = args.ped_file
#     discord_geno_file = args.discord_geno_file
#     hwe_file = args.hwe_file
#     gender_file = args.gender_info
#
#     global dp_threshold, gq_threshold
#     dp_threshold = args.dp
#     gq_threshold = args.gq
#
#     # discordant genotype file processing
#     discord_geno_dict = getDiscordInfo(discord_geno_file)
#
#     # vcf file processing and data recording
#     vcfProcessing(target_file, stat_file, gc_file, ped_file, discord_geno_dict, hwe_file, gender_file)