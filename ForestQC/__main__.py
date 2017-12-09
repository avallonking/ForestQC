import argparse, os
from ForestQC.stat import vcf_process, get_discord_info
from ForestQC.data_preprocessing import execute_split
from ForestQC.classification import execute_classification
from ForestQC.setOutlier import set_outlier
from ForestQC.get_ref_genome_GC import execute_compute_gc

def parse_args():
    # parse arguments
    parser = argparse.ArgumentParser(description='ForestQC: variant quality control based on random forest model')
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
    stat_parser.add_argument('-af', '--additional_features', default=None, type=str, required=False,
                             help='user-defined features in comma-separated file format [optional]')
    # arguments for splitting the variants data to good, bad and gray variants according to the statistics
    split_parser = subparsers.add_parser('split', help='variants splitting help')
    split_parser_required = split_parser.add_argument_group('required arguments')
    split_parser_required.add_argument('-i', '--input', dest='raw_stat_file', required=True,
                                       help='a raw variant dataset for splitting')
    split_parser.add_argument('-o', '--output', dest='output_handle', default=None, required=False,
                              help='the latter part of output filename [optional]')
    split_parser.add_argument('-m', '--model', dest='model', required=False, default='B',
                              help='splitting method for a specific random forest model, A or B [optional]. default: B')
    split_parser.add_argument('-af', '--additional_features', required=False, default=None, dest='user_feature_names',
                              help='names of additional user-defined features, can be a tab-separated text file or '
                                   'a comma-separated string')
    split_parser.add_argument('-t', '--thresholds', required=False, default=None, dest='thresholds_setting', type=str,
                              help='thresholds for each filter using in file splitting, should be a text file')

    # arguments for variant classification
    classify_parser = subparsers.add_parser('classify', help='variants classification help')
    classify_parser_required = classify_parser.add_argument_group('required arguments')
    classify_parser_required.add_argument('-g', '--good', dest='good_var', required=True, help='good variants')
    classify_parser_required.add_argument('-b', '--bad', dest='bad_var', required=True, help='bad variants')
    classify_parser_required.add_argument('-y', '--gray', dest='gray_var', required=True, help='gray variants')
    classify_parser.add_argument('-m', '--model', dest='model', required=False, default='B',
                                 help='type of random forest model, A or B [optional]. default: B')
    classify_parser.add_argument('-o', '--output', dest='output_suffix', required=False, default='variants.tsv',
                                 help='the latter part of output filename [optional]. default: variants.tsv')
    classify_parser.add_argument('-af', '--additional_features', dest='user_features', default=None, required=False,
                                 help='names of additional user-defined features used in the model, '
                                      'can be a tab-separated text file or a comma-separated string')

    # arguments for set outlier_gq and outlier_dp
    set_outlier_parser = subparsers.add_parser('set_outlier', help='set outlier help')
    set_outlier_parser_required = set_outlier_parser.add_argument_group('required arguments')
    set_outlier_parser_required.add_argument('-i', '--input', required=True, dest='file_list',
                                            help='input file(s). please separate the files with comma if there'
                                                 'are multiple files')

    # arguments for calculate GC content for reference genome
    compute_gc_parser = subparsers.add_parser('compute_gc', help='compute GC content help')
    compute_gc_parser_required = compute_gc_parser.add_argument_group('required arguments')
    compute_gc_parser_required.add_argument('-i', '--input', required=True, dest='ref_genome',
                                            help='reference genome fasta file')
    compute_gc_parser_required.add_argument('-o', '--output', required=True, dest='out_gc',
                                            help='output filename')
    compute_gc_parser.add_argument('-s', '--window_size', dest='window_size', type=int, default=1000,
                                   help='sliding window size for calculating GC content, default: 1000')


    args = vars(parser.parse_args())
    command = args.pop('command', None)
    return command, args

def main_set_outlier(**kwargs):
    file_list = kwargs['file_list'].split(',')
    set_outlier(file_list)

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

    # discordant genotype file processing
    discord_geno_dict = get_discord_info(discord_geno_file)
    # vcf file processing and data recording
    vcf_process(target_file, stat_file, gc_file, ped_file, discord_geno_dict, hwe_file, gender_file, dp, gq)

def main_split(**kwargs):
    input_file = kwargs['raw_stat_file']
    output_file = kwargs['output_handle']
    model = kwargs['model']
    user_feature_names_unprocessed = kwargs['user_feature_names']
    thresholds_setting = kwargs['thresholds_setting']

    if not user_feature_names_unprocessed:
        if os.path.isfile(user_feature_names_unprocessed):
            with open(user_feature_names_unprocessed, 'r') as u:
                line = u.readline()
                user_feature_names = line.strip().split('\t')
        else:
            user_feature_names = user_feature_names_unprocessed.strip().split(',')
    else:
        user_feature_names = None

    execute_split(input_file, output_file, model, user_feature_names, thresholds_setting)

def main_classify(**kwargs):
    good_var = kwargs['good_var']
    bad_var = kwargs['bad_var']
    gray_var = kwargs['gray_var']
    model = kwargs['model']
    output_suffix = kwargs['output_suffix']
    user_features_unprocessed = kwargs['user_features']

    if not user_features_unprocessed:
        if os.path.isfile(user_features_unprocessed):
            with open(user_features_unprocessed, 'r') as u:
                line = u.readline()
                user_features = line.strip().split('\t')
        else:
            user_features = user_features_unprocessed.strip().split(',')
    else:
        user_features = None

    execute_classification(good_var, bad_var, gray_var, model, output_suffix, user_features)

def main_compute_gc(**kwargs):
    ref = kwargs['ref_genome']
    out = kwargs['out_gc']
    window_size = kwargs['window_size']
    execute_compute_gc(ref, out, window_size)

def main():
    command_functions = {'stat': main_stat, 'split': main_split, 'classify': main_classify,
                         'set_outlier': main_set_outlier, 'compute_gc': main_compute_gc}
    command, args = parse_args()
    command_functions[command](**args)

if __name__ == '__main__':
    main()