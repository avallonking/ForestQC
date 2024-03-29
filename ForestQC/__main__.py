import argparse
import os
from ForestQC.stat import vcf_process, get_discord_info
from ForestQC.data_preprocessing import execute_split
from ForestQC.classification import execute_classification
from ForestQC.setOutlier import set_outlier
from ForestQC.get_ref_genome_GC import execute_compute_gc

def parse_args():
    # parse arguments
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='ForestQC command help', dest='command')
    subparsers.required = True

    # arguments for calculating statistics for variants
    stat_parser = subparsers.add_parser('stat', help='variants statistics help')
    stat_parser_required = stat_parser.add_argument_group('required arguments')
    stat_parser_required.add_argument('-i', '--input', dest='target_file', required=True,
                                      help='input vcf or vcf.gz file')
    stat_parser_required.add_argument('-o', '--output', dest='stat_file', required=True, help='output file name')
    stat_parser_required.add_argument('-c', '--gc', dest='gc_file', required=True,
                                      help='GC content of reference genome')
    stat_parser.add_argument('-p', '--ped', dest='ped_file', required=False, default=None,
                             help='pedigree file [optional]')
    stat_parser.add_argument('-d', '--discord_geno_file', required=False, dest='discord_geno_file', default=None,
                             help='discordant genotype file [optional]')
    stat_parser.add_argument('-w', '--hwe', dest='hwe_file', required=False, default=None,
                             help='HWE p-value file [optional]')
    stat_parser.add_argument('-g', '--gender', dest='gender_info', required=False, default=None,
                             help='information about gender of each individual [optional]')
    stat_parser.add_argument('--dp', dest='dp', default=34, type=float, required=False,
                             help='depth threshold [optional], default = 34')
    stat_parser.add_argument('--gq', dest='gq', default=99, type=float, required=False,
                             help='genotype quality threshold [optional], default = 99')
    stat_parser.add_argument('-as', '--additional_statistics', default=None, required=False, dest='add_stat_file',
                             help='user-defined features or filters in tab-separated file format [optional]')

    # arguments for splitting the variants data to good, bad and gray variants according to the statistics
    split_parser = subparsers.add_parser('split', help='variants splitting help')
    split_parser_required = split_parser.add_argument_group('required arguments')
    split_parser_required.add_argument('-i', '--input', dest='raw_stat_file', required=True,
                                       help='a raw variant dataset for splitting')
    split_parser.add_argument('-o', '--output', dest='output_handle', default=None, required=False,
                              help='the latter part of output filename [optional]')
    split_parser.add_argument('-m', '--model', dest='model', required=False, default='B',
                              help='splitting method for a specific random forest model, A or B [optional]. default: B')
    split_parser.add_argument('-as', '--additional_statistics', required=False, default=None, dest='user_stat_names',
                              help='names of additional user-defined features or filters added in the \"stat\" step, can'
                                   ' be a tab-separated text file (the file in the \"stat\" step) or a comma-separated '
                                   'string [optional]')
    split_parser.add_argument('-t', '--thresholds', required=False, default=None, dest='thresholds_setting', type=str,
                              help='thresholds for each filter using in file splitting, should be a tab-separated '
                                   'file. if user-defined features are included in statistics files, this option is '
                                   'required.')
    # split_parser.add_argument('-f', '--filters', dest='filters', default=None, required=False,
    #                              help='names of filters to be used. it should be a comma-separated string [optional]. '
    #                                   'default: Mendel_Error,Missing_Rate,HWE,ABHet')

    # arguments for variant classification
    classify_parser = subparsers.add_parser('classify', help='variants classification help')
    classify_parser_required = classify_parser.add_argument_group('required arguments')
    classify_parser_required.add_argument('-g', '--good', dest='good_var', required=True, help='good variants')
    classify_parser_required.add_argument('-b', '--bad', dest='bad_var', required=True, help='bad variants')
    classify_parser_required.add_argument('-y', '--gray', dest='gray_var', required=True, help='gray variants')
    classify_parser.add_argument('-t', '--threshold', dest='prob_threshold', required=False, default=0.50, type = float,
                                 help='probability threshold for variant classification, default = 0.50')
    classify_parser.add_argument('-m', '--model', dest='model', required=False, default='B',
                                 help='type of random forest model, A or B [optional]. default: B')
    classify_parser.add_argument('-o', '--output', dest='output_suffix', required=False, default='variants.tsv',
                                 help='the latter part of output filename [optional]. default: variants.tsv')
    # classify_parser.add_argument('-af', '--additional_features', dest='user_features', default=None, required=False,
    #                              help='names of additional user-defined features used in the model, '
    #                                   'can be a tab-separated text file or a comma-separated string.')
    classify_parser.add_argument('-f', '--features', dest='features', default='Mean_DP,Mean_GQ,SD_DP,SD_GQ,Outlier_DP'
                                                                              ',Outlier_GQ,GC', required=False,
                                 help='names of features to be used in the random forest model. it should be a '
                                      'comma-separated string [optional]. '
                                      'default: Mean_DP,Mean_GQ,SD_DP,SD_GQ,Outlier_DP,Outlier_GQ,GC')

    # arguments for set outlier_gq and outlier_dp
    set_outlier_parser = subparsers.add_parser('set_outlier', help='set outlier help')
    set_outlier_parser_required = set_outlier_parser.add_argument_group('required arguments')
    set_outlier_parser_required.add_argument('-i', '--input', required=True, dest='file_list',
                                            help='input file(s). please separate the filenames with comma if there'
                                                 'are multiple files')
    set_outlier_parser.add_argument('-m', '--mem', required=False, dest='mem', default='2g',
                                    help='memory usage for external sort. default: 2G')
    set_outlier_parser.add_argument('-t', '--temp', required=False, dest='temp_dir', default=None,
                                    help='directory to store intermediate files. default: current working directory')

    # arguments for calculate GC content for reference genome
    compute_gc_parser = subparsers.add_parser('compute_gc', help='compute GC content help')
    compute_gc_parser_required = compute_gc_parser.add_argument_group('required arguments')
    compute_gc_parser_required.add_argument('-i', '--input', required=True, dest='ref_genome',
                                            help='reference genome fasta file')
    compute_gc_parser_required.add_argument('-o', '--output', required=True, dest='out_gc',
                                            help='output filename')
    compute_gc_parser.add_argument('-s', '--window_size', dest='window_size', type=int, default=1000, required=False,
                                   help='sliding window size for calculating GC content [optional], default: 1000')


    args = vars(parser.parse_args())
    command = args.pop('command', None)
    return command, args

def main_set_outlier(**kwargs):
    file_list = kwargs['file_list'].split(',')
    mem = kwargs['mem']
    temp_dir = kwargs['temp_dir']
    if temp_dir == None:
        temp_dir = os.getcwd()
    set_outlier(file_list, temp_dir, 'temp.external_sort.out', mem)

def main_stat(**kwargs):
    target_file = kwargs['target_file']
    stat_file = kwargs['stat_file']
    gc_file = kwargs['gc_file']
    ped_file = kwargs['ped_file']
    discord_geno_file = kwargs['discord_geno_file']
    hwe_file = kwargs['hwe_file']
    gender_file = kwargs['gender_info']
    dp = kwargs['dp']
    gq = kwargs['gq']
    feature_file = kwargs['add_stat_file']

    # discordant genotype file processing
    discord_geno_dict = get_discord_info(discord_geno_file)
    # vcf file processing and data recording
    vcf_process(target_file, stat_file, gc_file, ped_file, discord_geno_dict, hwe_file, gender_file, dp, gq,
                feature_file)

def main_split(**kwargs):
    input_file = kwargs['raw_stat_file']
    output_file = kwargs['output_handle']
    model = kwargs['model']
    user_feature_names_unprocessed = kwargs['user_stat_names']
    thresholds_setting = kwargs['thresholds_setting']
    # filters = kwargs['filters']
    #
    # try:
    #     filters = filters.strip().split(',')
    # except TypeError:
    #     filters = None

    try:
        if os.path.isfile(user_feature_names_unprocessed):
            with open(user_feature_names_unprocessed, 'r') as u:
                line = u.readline()
                user_feature_names = line.strip().split('\t')
                try:
                    user_feature_names.remove('RSID')
                except ValueError:
                    pass
        else:
            user_feature_names = user_feature_names_unprocessed.strip().split(',')
    except TypeError:
        user_feature_names = None

    execute_split(input_file, output_file, model, user_feature_names, thresholds_setting)

def main_classify(**kwargs):
    good_var = kwargs['good_var']
    bad_var = kwargs['bad_var']
    gray_var = kwargs['gray_var']
    model = kwargs['model']
    output_suffix = kwargs['output_suffix']
    features = kwargs['features']
    prob_threshold = kwargs['prob_threshold']

    features = features.strip().split(',')
    # user_features_unprocessed = kwargs['user_features']
    # try:
    #     if os.path.isfile(user_features_unprocessed):
    #         with open(user_features_unprocessed, 'r') as u:
    #             line = u.readline()
    #             user_features = line.strip().split('\t')
    #             try:
    #                 user_features.remove('RSID')
    #             except ValueError:
    #                 pass
    #     else:
    #         user_features = user_features_unprocessed.strip().split(',')
    # except TypeError:
    #     user_features = None

    execute_classification(good_var, bad_var, gray_var, model, output_suffix, features, prob_threshold)

def main_compute_gc(**kwargs):
    ref = kwargs['ref_genome']
    out = kwargs['out_gc']
    window_size = kwargs['window_size']
    execute_compute_gc(ref, out, window_size)

def main():
    print('ForestQC v1.1.5.7 by Jae Hoon Sul Lab at UCLA')
    print('--Quality control on genetic variants from next-generation sequencing data using random forest')
    print()
    command_functions = {'stat': main_stat, 'split': main_split, 'classify': main_classify,
                         'set_outlier': main_set_outlier, 'compute_gc': main_compute_gc}
    command, args = parse_args()
    command_functions[command](**args)

if __name__ == '__main__':
    main()
