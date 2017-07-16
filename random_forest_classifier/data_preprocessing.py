# This program is for the preparation of traning set and test set. It divides the variants
# dataset into 3 datasets containing good variants, bad variants and grey variants, respectively.

import pandas as pd
import sys

def separateDataA(variants):
    # input: all variants with statistics, a pandas dataframe
    # output: good, bad and grey variants

    # columns names
    columns = ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom']
    variants.columns = columns

    # good variants
    good = variants[~(variants['Mendel_Error'] > 0)]
    good = good[good['Missing_Rate'] < 0.005]
    good = good[good['HWE'] > 0.01]
    good = good[~(good['Discordant_Geno'] > 1)]
    good['Good'] = 1

    # bad variants
    rare_variants = variants[variants['MAF'] < 0.03]
    common_variants = variants[variants['MAF'] >= 0.03]
    rare_bad = rare_variants[~(rare_variants['Mendel_Error'] <= 3)]
    rare_bad = rare_bad[~(rare_bad['Discordant_Geno'] <= 6)]
    rare_bad = rare_bad[rare_bad['Missing_Rate'] > 0.02]
    rare_bad = rare_bad[rare_bad['HWE'] < 5e-3]
    common_bad = common_variants[~(common_variants['Mendel_Error'] <= 5)]
    common_bad = common_bad[~(common_bad['Discordant_Geno'] <= 6)]
    common_bad = common_bad[common_bad['Missing_Rate'] > 0.03]
    common_bad = common_bad[common_bad['HWE'] < 5e-4]

    bad = pd.concat([rare_bad, common_bad])
    bad.drop_duplicates(inplace=True)
    bad['Good'] = 0

    # grey variants
    grey = variants[~variants['RSID'].isin(good['RSID'])]
    grey = grey[~grey['RSID'].isin(bad['RSID'])]

    return good, bad, grey

def separateDataB(variants):
    # input: all variants with statistics, a pandas dataframe
    # output: good, bad and grey variants

    # columns names
    columns = ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom']
    variants.columns = columns

    # good variants
    good = variants[~(variants['Mendel_Error'] > 0)]
    good = good[good['Missing_Rate'] < 0.005]
    good = good[good['HWE'] > 0.01]
    good = good[~(good['Discordant_Geno'] > 1)]
    good = good[~good['ABHet'].isnull()]
    good['Good'] = 1

    # bad variants
    rare_variants = variants[variants['MAF'] < 0.03]
    common_variants = variants[variants['MAF'] >= 0.03]
    rare_bad = rare_variants[~(rare_variants['Mendel_Error'] <= 3)]
    rare_bad = rare_bad[~(rare_bad['Discordant_Geno'] <= 6)]
    rare_bad = rare_bad[rare_bad['Missing_Rate'] > 0.02]
    rare_bad = rare_bad[rare_bad['HWE'] < 5e-3]
    common_bad = common_variants[~(common_variants['Mendel_Error'] <= 5)]
    common_bad = common_bad[~(common_bad['Discordant_Geno'] <= 6)]
    common_bad = common_bad[common_bad['Missing_Rate'] > 0.03]
    common_bad = common_bad[common_bad['HWE'] < 5e-4]

    outliers_1 = variants[variants['Discordant_Geno'] > 20]
    outliers_2 = pd.concat([rare_variants[rare_variants['Mendel_Error'] > 8], common_variants[common_variants['Mendel_Error'] > 10]])
    outliers_3 = pd.concat([rare_variants[rare_variants['Missing_Rate'] > 0.08], common_variants[common_variants['Missing_Rate'] > 0.10]])
    outliers_4 = pd.concat([rare_variants[rare_variants['HWE'] < 1e-3], common_variants[common_variants['HWE'] < 1e-8]])

    bad = pd.concat([rare_bad, common_bad, outliers_1, outliers_2, outliers_3, outliers_4])
    bad.drop_duplicates(inplace=True)
    bad = bad[~bad['ABHet'].isnull()]
    bad['Good'] = 0

    # grey variants
    grey = variants[~variants['RSID'].isin(good['RSID'])]
    grey = grey[~grey['RSID'].isin(bad['RSID'])]
    grey = grey[~grey['ABHet'].isnull()]

    return good, bad, grey

def separateDataC(variants):
    # input: all variants with statistics, a pandas dataframe
    # output: good, bad and grey variants

    # columns names
    columns = ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom']
    variants.columns = columns

    # good variants
    good = variants[~(variants['Mendel_Error'] > 0)]
    good = good[good['Missing_Rate'] < 0.005]
    good = good[good['HWE'] > 0.01]
    good = good[abs(good['ABHet'] - 0.5) <= 0.2]
    good = good[~(good['Discordant_Geno'] > 1)]
    good['Good'] = 1

    # bad variants
    rare_variants = variants[variants['MAF'] < 0.03]
    common_variants = variants[variants['MAF'] >= 0.03]
    rare_bad = rare_variants[~(rare_variants['Mendel_Error'] <= 3)]
    rare_bad = rare_bad[~(rare_bad['Discordant_Geno'] <= 6)]
    rare_bad = rare_bad[rare_bad['Missing_Rate'] > 0.02]
    rare_bad = rare_bad[rare_bad['HWE'] < 5e-3]
    rare_bad = rare_bad[abs(rare_bad['ABHet'] - 0.5) >= 0.25]
    common_bad = common_variants[~(common_variants['Mendel_Error'] <= 5)]
    common_bad = common_bad[~(common_bad['Discordant_Geno'] <= 6)]
    common_bad = common_bad[common_bad['Missing_Rate'] > 0.03]
    common_bad = common_bad[common_bad['HWE'] < 5e-4]
    common_bad = common_bad[abs(common_bad['ABHet'] - 0.5) >= 0.25]

    outliers_1 = variants[variants['Discordant_Geno'] > 20]
    outliers_2 = pd.concat([rare_variants[rare_variants['Mendel_Error'] > 8], common_variants[common_variants['Mendel_Error'] > 10]])
    outliers_3 = pd.concat([rare_variants[rare_variants['Missing_Rate'] > 0.08], common_variants[common_variants['Missing_Rate'] > 0.10]])
    outliers_4 = pd.concat([rare_variants[rare_variants['HWE'] < 1e-3], common_variants[common_variants['HWE'] < 1e-8]])
    # outliers_5 = pd.concat([rare_variants[abs(rare_variants['ABHet'] - 0.5) >= 0.25], common_variants[abs(common_variants['ABHet'] - 0.5) >= 0.25]])

    # bad = pd.concat([rare_bad, common_bad, outliers_1, outliers_2, outliers_3, outliers_4, outliers_5])
    bad = pd.concat([rare_bad, common_bad, outliers_1, outliers_2, outliers_3, outliers_4])
    bad.drop_duplicates(inplace=True)
    bad['Good'] = 0

    # grey variants
    grey = variants[~variants['RSID'].isin(good['RSID'])]
    grey = grey[~grey['RSID'].isin(bad['RSID'])]

    return good, bad, grey

def preprocessing(data):
  columns = ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom']
  data.columns = columns
  data = data[~data['Mean_DP'].isnull()]
  data = data[~data['Mean_GQ'].isnull()]
  data = data[~data['SD_DP'].isnull()]
  data = data[~data['SD_GQ'].isnull()]
  data = data[~data['Outlier_DP'].isnull()]
  data = data[~data['Outlier_GQ'].isnull()]
  return data

def main():
  input_file = sys.argv[1]
  output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.split('/')[-1]
  try:
    data = pd.read_table(input_file, header=None)
  except pd.io.common.EmptyDataError:
    exit(1)
  data = preprocessing(data)
  good, bad, grey = separateDataB(data)
  good.to_csv('good.' + output_file, index=False, header=False, sep='\t', na_rep='NA')
  bad.to_csv('bad.' + output_file, index=False, header=False, sep='\t', na_rep='NA')
  grey.to_csv('grey.' + output_file, index=False, header=False, sep='\t', na_rep='NA')

if __name__ == '__main__':
  main()
