# This program is for the preparation of traning set and test set. It divides the variants
# dataset into 3 datasets containing good variants, bad variants and grey variants, respectively.

import pandas as pd
import os

# deprecated
# def separateDataA(variants):
#     # input: all variants with statistics, a pandas dataframe
#     # output: good, bad and grey variants
#
#     # columns names
#     columns = ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom', 'GC']
#     variants.columns = columns
#
#     # good variants
#     good = variants[~(variants['Mendel_Error'] > 0)]
#     good = good[good['Missing_Rate'] < 0.005]
#     good = good[good['HWE'] > 0.01]
#     good = good[~(good['Discordant_Geno'] > 1)]
#     good['Good'] = 1
#
#     # bad variants
#     rare_variants = variants[variants['MAF'] < 0.03]
#     common_variants = variants[variants['MAF'] >= 0.03]
#     rare_bad = rare_variants[~(rare_variants['Mendel_Error'] <= 3)]
#     rare_bad = rare_bad[~(rare_bad['Discordant_Geno'] <= 6)]
#     rare_bad = rare_bad[rare_bad['Missing_Rate'] > 0.02]
#     rare_bad = rare_bad[rare_bad['HWE'] < 5e-3]
#     common_bad = common_variants[~(common_variants['Mendel_Error'] <= 5)]
#     common_bad = common_bad[~(common_bad['Discordant_Geno'] <= 6)]
#     common_bad = common_bad[common_bad['Missing_Rate'] > 0.03]
#     common_bad = common_bad[common_bad['HWE'] < 5e-4]
#
#     bad = pd.concat([rare_bad, common_bad])
#     bad.drop_duplicates(inplace=True)
#     bad['Good'] = 0
#
#     # grey variants
#     grey = variants[~variants['RSID'].isin(good['RSID'])]
#     grey = grey[~grey['RSID'].isin(bad['RSID'])]
#
#     return good, bad, grey

def separateDataA(variants, thresholds_setting):
    # input: all variants with statistics, a pandas dataframe
    # output: good, bad and grey variants

    variants.loc[variants['ABHet'].isnull(), 'ABHet'] = variants['ABHet'].median()
    variants.loc[variants['ABHom'].isnull(), 'ABHom'] = variants['ABHom'].median()

    good_thresholds, bad_thresholds, outlier_thresholds = set_thresholds(thresholds_setting)

    # good variants
    good = variants[~(variants['Mendel_Error'] > good_thresholds['Mendel_Error'])]
    good = good[good['Missing_Rate'] < good_thresholds['Missing_Rate']]
    good = good[good['HWE'] > good_thresholds['HWE']]
    # good = good[~(good['Discordant_Geno'] > good_thresholds['Discordant_Geno'])]
    good['Good'] = 1

    # bad variants
    rare_variants = variants[variants['MAF'] < bad_thresholds['all']['MAF']]
    common_variants = variants[variants['MAF'] >= bad_thresholds['all']['MAF']]
    rare_bad = rare_variants[~(rare_variants['Mendel_Error'] <= bad_thresholds['rare']['Mendel_Error'])]
    # rare_bad = rare_bad[~(rare_bad['Discordant_Geno'] <= bad_thresholds['rare']['Discordant_Geno'])]
    rare_bad = rare_bad[rare_bad['Missing_Rate'] > bad_thresholds['rare']['Missing_Rate']]
    rare_bad = rare_bad[rare_bad['HWE'] < bad_thresholds['rare']['HWE']]
    common_bad = common_variants[~(common_variants['Mendel_Error'] <= bad_thresholds['common']['Mendel_Error'])]
    # common_bad = common_bad[~(common_bad['Discordant_Geno'] <= bad_thresholds['common']['Discordant_Geno'])]
    common_bad = common_bad[common_bad['Missing_Rate'] > bad_thresholds['common']['Missing_Rate']]
    common_bad = common_bad[common_bad['HWE'] < bad_thresholds['common']['HWE']]

    # outliers_1 = pd.concat([rare_variants[rare_variants['Discordant_Geno'] > outlier_thresholds['rare']['Discordant_Geno']],
                            # common_variants[common_variants['Discordant_Geno'] > outlier_thresholds['common']['Discordant_Geno']]])
    outliers_2 = pd.concat([rare_variants[rare_variants['Mendel_Error'] > outlier_thresholds['rare']['Mendel_Error']],
                            common_variants[common_variants['Mendel_Error'] > outlier_thresholds['common']['Mendel_Error']]])
    outliers_3 = pd.concat([rare_variants[rare_variants['Missing_Rate'] > outlier_thresholds['rare']['Missing_Rate']],
                            common_variants[common_variants['Missing_Rate'] > outlier_thresholds['common']['Missing_Rate']]])
    outliers_4 = pd.concat([rare_variants[rare_variants['HWE'] < outlier_thresholds['rare']['HWE']],
                            common_variants[common_variants['HWE'] < outlier_thresholds['common']['HWE']]])

    bad = pd.concat([rare_bad, common_bad, outliers_2, outliers_3, outliers_4])
    bad.drop_duplicates(inplace=True)
    bad['Good'] = 0

    # grey variants
    grey = variants[~variants['RSID'].isin(good['RSID'])]
    grey = grey[~grey['RSID'].isin(bad['RSID'])]

    return good, bad, grey

def separateDataB(variants, thresholds_setting):
    # input: all variants with statistics, a pandas dataframe
    # output: good, bad and grey variants

    good_thresholds, bad_thresholds, outlier_thresholds = set_thresholds(thresholds_setting)

    # good variants
    good = variants[~(variants['Mendel_Error'] > good_thresholds['Mendel_Error'])]
    good = good[good['Missing_Rate'] < good_thresholds['Missing_Rate']]
    good = good[good['HWE'] > good_thresholds['HWE']]
    # good = good[~(good['Discordant_Geno'] > good_thresholds['Discordant_Geno'])]
    good = good[abs(good['ABHet'] - 0.5) <= good_thresholds['ABHet_deviation']]
    good['Good'] = 1

    # bad variants
    rare_variants = variants[variants['MAF'] < bad_thresholds['all']['MAF']]
    common_variants = variants[variants['MAF'] >= bad_thresholds['all']['MAF']]

    rare_bad = rare_variants[~(rare_variants['Mendel_Error'] <= bad_thresholds['rare']['Mendel_Error'])]
    # rare_bad = rare_bad[~(rare_bad['Discordant_Geno'] <= bad_thresholds['rare']['Discordant_Geno'])]
    rare_bad = rare_bad[rare_bad['Missing_Rate'] > bad_thresholds['rare']['Missing_Rate']]
    rare_bad = rare_bad[rare_bad['HWE'] < bad_thresholds['rare']['HWE']]
    rare_bad = rare_bad[abs(rare_bad['ABHet'] - 0.5) >= bad_thresholds['rare']['ABHet_deviation']]

    common_bad = common_variants[~(common_variants['Mendel_Error'] <= bad_thresholds['common']['Mendel_Error'])]
    # common_bad = common_bad[~(common_bad['Discordant_Geno'] <= bad_thresholds['common']['Discordant_Geno'])]
    common_bad = common_bad[common_bad['Missing_Rate'] > bad_thresholds['common']['Missing_Rate']]
    common_bad = common_bad[common_bad['HWE'] < bad_thresholds['common']['HWE']]
    common_bad = common_bad[abs(common_bad['ABHet'] - 0.5) >= bad_thresholds['common']['ABHet_deviation']]

    # outliers_1 = pd.concat([rare_variants[rare_variants['Discordant_Geno'] > outlier_thresholds['rare']['Discordant_Geno']],
                            # common_variants[common_variants['Discordant_Geno'] > outlier_thresholds['common']['Discordant_Geno']]])
    outliers_2 = pd.concat([rare_variants[rare_variants['Mendel_Error'] > outlier_thresholds['rare']['Mendel_Error']],
                            common_variants[common_variants['Mendel_Error'] > outlier_thresholds['common']['Mendel_Error']]])
    outliers_3 = pd.concat([rare_variants[rare_variants['Missing_Rate'] > outlier_thresholds['rare']['Missing_Rate']],
                            common_variants[common_variants['Missing_Rate'] > outlier_thresholds['common']['Missing_Rate']]])
    outliers_4 = pd.concat([rare_variants[rare_variants['HWE'] < outlier_thresholds['rare']['HWE']],
                            common_variants[common_variants['HWE'] < outlier_thresholds['common']['HWE']]])

    bad = pd.concat([rare_bad, common_bad, outliers_2, outliers_3, outliers_4])
    bad.drop_duplicates(inplace=True)
    bad['Good'] = 0

    # grey variants
    grey = variants[~variants['RSID'].isin(good['RSID'])]
    grey = grey[~grey['RSID'].isin(bad['RSID'])]

    return good, bad, grey

def preprocessing(data, user_feature_names):
    try:
        columns = ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP',
                   'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom', 'GC'] \
                  + user_feature_names
        assert len(columns) == data.shape[1], 'Missing names of user-defined features'
        data.columns = columns
        # impute missing values
        for col in user_feature_names:
            data.loc[data[col].isnull(), col] = data[col].median()
    except TypeError:
        columns = ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP',
                   'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom', 'GC']
        assert len(columns) == data.shape[1], 'Missing names of user-defined features'
        data.columns = columns

    data = data[~data['MAF'].isnull()]
    data.loc[data['Mean_DP'].isnull(), 'Mean_DP'] = data['Mean_DP'].median()
    data.loc[data['Mean_GQ'].isnull(), 'Mean_GQ'] = data['Mean_GQ'].median()
    data.loc[data['SD_DP'].isnull(), 'SD_DP'] = data['SD_DP'].median()
    data.loc[data['SD_GQ'].isnull(), 'SD_GQ'] = data['SD_GQ'].median()
    data.loc[data['Outlier_DP'].isnull(), 'Outlier_DP'] = data['Outlier_DP'].median()
    data.loc[data['Outlier_GQ'].isnull(), 'Outlier_GQ'] = data['Outlier_GQ'].median()
    # data.loc[data['ABHom'].isnull(), 'ABHom'] = data['ABHom'].median()
    # data.loc[data['ABHet'].isnull(), 'ABHet'] = data['ABHet'].median()
    data.loc[data['GC'].isnull(), 'GC'] = data['GC'].median()
    return data

def print_thresholds(thresholds, type):
    # thresholds is a one-demensional dict
    for _filter, _threshold in thresholds.items():
        if type == 'good':
            if _filter == 'ABHet_deviation':
                print('{} <= {} <= {}'.format(0.5 - _threshold, _filter, 0.5 + _threshold))
            elif _filter == 'Missing_Rate':
                print('{} < {}'.format(_filter, _threshold))
            elif _filter == 'HWE':
                print('{} > {}'.format(_filter, _threshold))
            else:
                print('{} <= {}'.format(_filter, _threshold))
        else:
            if _filter == 'ABHet_deviaiton':
                print('\t{} >= {} or {} <= {}'.format(_filter, 0.5 + _threshold, _filter, 0.5 - _threshold))
            elif _filter == 'HWE':
                print('\t{} < {}'.format(_filter, _threshold))
            else:
                print('\t{} > {}'.format(_filter, _threshold))

def set_thresholds(thresholds_setting):
    # input format (tab-separated):
    #   good    all filter  threshold
    #   bad all MAF threshold
    #   bad rare    Mendel_Error   threshold
    #   bad common  Missing_Rate    threshold
    #   outlier rare    Missing_Rate   threshold
    #   outlier common  Missing_Rate    threshold

    good_thresholds = {'Mendel_Error': round(3 / 446, 5), 'Missing_Rate':0.005, 'HWE': 0.01,
                       'ABHet_deviation': 0.20}
    bad_thresholds = {'all': {'MAF': 0.03},
                      'rare': {'Mendel_Error': round(3 / 446, 5), 'Missing_Rate': 0.02, 'HWE': 5e-3,
                      'ABHet_deviation': 0.25},
                      'common': {'Mendel_Error': round(5 / 446, 5), 'Missing_Rate': 0.03, 'HWE': 5e-4,
                      'ABHet_deviation': 0.25}
                     }
    outlier_thresholds = {'rare': {'Mendel_Error': round(8 / 446, 5), 'Missing_Rate': 0.08, 'HWE': 1e-3},
                          'common': {'Mendel_Error': round(10 / 446, 5), 'Missing_Rate': 0.10, 'HWE': 1e-8}
                          }
    try:
        with open(thresholds_setting, 'r') as t:
            for line in t:
                _info = line.strip().split('\t')
                _type = _info[0]
                _scope = _info[1]
                _filter = _info[2]
                threshold = _info[3]
                if _type == 'bad':
                    bad_thresholds[_scope][_filter] = float(threshold)
                elif _type == 'good':
                    good_thresholds[_filter] = float(threshold)
                elif _type == 'outlier':
                    outlier_thresholds[_scope][_filter] = float(threshold)
    except TypeError:
        pass
    finally:
        # show the thresholds setting
        print('\nCurrent filter settings:')
        print('\nGood variants')
        print('----------------')
        print_thresholds(good_thresholds, 'good')
        print('\nBad variants')
        print('----------------')
        print('Rare variants (MAF < {}):'.format(bad_thresholds['all']['MAF']))
        print_thresholds(bad_thresholds['rare'], 'bad')
        print('\nCommon variants (MAF >= {}):'.format(bad_thresholds['all']['MAF']))
        print_thresholds(bad_thresholds['common'], 'bad')
        print('\nOutlier variants')
        print('----------------')
        print('Rare variants (MAF < {}):'.format(bad_thresholds['all']['MAF']))
        print_thresholds(outlier_thresholds['rare'], 'outlier')
        print('\nCommon variants (MAF >= {}):'.format(bad_thresholds['all']['MAF']))
        print_thresholds(outlier_thresholds['common'], 'outlier')

        return good_thresholds, bad_thresholds, outlier_thresholds

def execute_split(input_file, output_file, model, user_feature_names, thresholds_setting):
    # user_feature_names should be a list
    if not output_file:
        output_file = input_file.split('/')[-1]
    try:
        print('Loading data...')
        data = pd.read_table(input_file, header=None)
    except pd.io.common.EmptyDataError:
        raise FileNotFoundError('No such file: ' + input_file + '\n')
        exit(1)

    print('Data processing...')
    data = preprocessing(data, user_feature_names)

    model_selection = {'A': separateDataA, 'B': separateDataB}
    good, bad, grey = model_selection[model](data, thresholds_setting)

    print('\nNumber of variants')
    print('Good variants: ' + str(good.shape[0]))
    print('Bad variants: ' + str(bad.shape[0]))
    print('Grey variants: ' + str(grey.shape[0]))

    print('\nWriting data...')
    dir = os.path.dirname(os.path.realpath(input_file))
    good.to_csv(dir + '/' + 'good.' + output_file, index=False, sep='\t', na_rep='NA')
    bad.to_csv(dir + '/' + 'bad.' + output_file, index=False, sep='\t', na_rep='NA')
    grey.to_csv(dir + '/' + 'gray.' + output_file, index=False, sep='\t', na_rep='NA')
    print('Done.')
