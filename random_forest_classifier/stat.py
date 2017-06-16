import gzip
import argparse
from vcf_stat import *

# Deprecated
# import sys
# from data_preprocessing import *
# from classification import *

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
    except OSError:
      pass
    return discord_geno_dict

def vcfProcessing(vcf_file, stat_file, ped_file, discord_geno_dict, hwe_file):
    f = gzip.open(vcf_file, 'rt') if 'gz' in vcf_file else open(vcf_file, 'r')
    o = open(stat_file, 'w')
    for line in f:
        if line.startswith('#CHROM'):
          sample_list = line.strip().split('\t')[DP_GQ_START_IDX:]
          break
    
    relationship = getFamilyRelation(ped_file, sample_list)
    control_samples_idx = getControlSamples(ped_file, sample_list)
    hwe_info = getHWE_Direct(hwe_file)
    
    for line2 in f:
        if line2.startswith('chr'):
          site_info = line2.split('\t')
          chr = site_info[0]
          pos = site_info[1]
          ref = site_info[3]
          alt = site_info[4]
          rsid = chr + ':' + pos
          maf = getMAF(line2)
          mean_dp, mean_gq, sd_dp, sd_gq, outlier_dp, outlier_gq = statDPGQ(line2)
          discordant_geno = getDiscordantGenotype(line2, discord_geno_dict)
          mendel_error = getMendel(line2, sample_list, relationship)
          missing_rate = getMissing(line2)
          hwe = getHWE(line2, control_samples_idx) if not hwe_info else hwe_info[rsid]
          abhet, abhom = getAB(line2)
          o.write('\t'.join([rsid, chr, pos, ref, alt, str(maf), str(mean_dp), str(mean_gq), str(sd_dp), str(sd_gq), str(outlier_dp), str(outlier_gq), str(discordant_geno), str(mendel_error), str(missing_rate), str(hwe), str(abhet), str(abhom)]) + '\n')
    o.close()
    f.close()

def main():
    # target_file: the vcf file for processing
    # stat_file: the output filename
    # ped_file: the pedigree file
    # discord_geno_file: the file that has the information about discordant genotype
    parser = argparse.ArgumentParser(description='Calculate statistics for one vcf file')
    parser.add_argument('-i', '--input', dest='target_file', help='Input vcf or vcf.gz file')
    parser.add_argument('-o', '--output', dest='stat_file', help='Output file name')
    parser.add_argument('-p', '--ped', dest='ped_file', default='NA', help='Pedigree file [optional]')
    parser.add_argument('-d', '--discord_geno_file', dest='discord_geno_file', default='NA', help='Discordant genotype file [optional]')
    parser.add_argument('-w', '--hwe', dest='hwe_file', default='NA', help='HWE p-value file [optional]')
    args = parser.parse_args()

    target_file = args.target_file
    stat_file = args.stat_file
    ped_file = args.ped_file
    discord_geno_file = args.discord_geno_file
    hwe_file = args.hwe_file

    # discordant genotype file processing
    discord_geno_dict = getDiscordInfo(discord_geno_file)

    # vcf file processing and data recording
    vcfProcessing(target_file, stat_file, ped_file, discord_geno_dict, hwe_file)

if __name__ == '__main__':
    main()
