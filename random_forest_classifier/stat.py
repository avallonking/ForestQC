import sys
import gzip
from vcf_stat import *

# Deprecated
# from data_preprocessing import *
# from classification import *

def getDiscordInfo(discord_geno_file):
    # discordant genotype file processing
    discord_geno_dict = {}
    d = open(discord_geno_file, 'r')
    for entry in d:
        snp_id = entry.split('\t')[0]
        discord_geno_num = entry.strip().split('\t')[1]
        discord_geno_dict[snp_id] = discord_geno_num
    d.close()
    return discord_geno_dict

def vcfProcessing(vcf_file, stat_file, ped_file, discord_geno_dict):
    f = gzip.open(vcf_file, 'rt') if 'gz' in vcf_file else open(vcf_file, 'r')
    o = open(stat_file, 'w')
    for line in f:
        if line.startswith('#CHROM'):
          sample_list = line.strip().split('\t')[DP_GQ_START_IDX:]
          break
    relationship = getFamilyRelation(ped_file, sample_list)
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
          discordant_geno = getDiscordantGenotype(line2, discord_geno_dict) if discord_geno_dict != 'NA' else 'NA'
          mendel_error = getMendel(line2, sample_list, relationship)
          missing_rate = getMissing(line2)
          hwe = getHWE(line2)
          abhet, abhom = getAB(line2)
          o.write('\t'.join([rsid, chr, pos, ref, alt, str(maf), str(mean_dp), str(mean_gq), str(sd_dp), str(sd_gq), str(outlier_dp), str(outlier_gq), str(discordant_geno), str(mendel_error), str(missing_rate), str(hwe), str(abhet), str(abhom)]) + '\n')
    o.close()
    f.close()

def main():
    # target_file: the vcf file for processing
    # stat_file: the output filename
    # ped_file: the pedigree file
    # discord_geno_file: the file that has the information about discordant genotype
    target_file = sys.argv[1]
    stat_file = sys.argv[2]
    ped_file = sys.argv[3]
    discord_geno_file = sys.argv[4] if len(sys.argv) >= 5 else 'NA'

    # discordant genotype file processing
    discord_geno_dict = getDiscordInfo(discord_geno_file) if discord_geno_file != 'NA' else 'NA'

    # vcf file processing and data recording
    vcfProcessing(target_file, stat_file, ped_file, discord_geno_dict)

if __name__ == '__main__':
    main()
