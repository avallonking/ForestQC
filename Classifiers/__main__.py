import argparse
from stat import *
from vcf_stat import *
from data_preprocessing import *
from classification import *
from setOutlier import *

def parse_args():
    # parse arguments
    parser = argparse.ArgumentParser(description='ForestQC: variant quality control based on random forest model')
    subparsers = parser.add_subparsers(help='ForestQC command help', dest='command')
    subparsers.required = True

    # arguments for calculating statistics for variants
    stat_parser = subparsers.add_parser('stat', help='variants statistics help')
    stat_parser_required = stat_parser.add_argument_group('required arguments')
    stat_parser_required.add_argument('-i', '--input', dest='target_file', required=True, help='input vcf or vcf.gz file')
    stat_parser_required.add_argument('-o', '--output', dest='stat_file', required=True, help='output file name')
    stat_parser_required.add_argument('-c', '--gc', dest='gc_file', required=True, help='GC content of reference genome')
    stat_parser.add_argument('-p', '--ped', dest='ped_file', required=False, default=None, help='pedigree file [optional]')
    stat_parser.add_argument('-d', '--discord_geno_file', required=False, dest='discord_geno_file', default=None,
                             help='discordant genotype file [optional]')
    stat_parser.add_argument('-w', '--hwe', dest='hwe_file', required=False, default=None, help='HWE p-value file [optional]')
    stat_parser.add_argument('-g', '--gender', dest='gender_info', required=False, default=None,
                             help='information about gender of each individual [optional]')
    stat_parser.add_argument('--dp', dest='dp', default=34, type=float, required=False, help='depth threshold [optional], default = 34')
    stat_parser.add_argument('--gq', dest='gq', default=99, type=float, required=False, help='genotype quality threshold [optional], default = 99')

    # arguments for splitting the variants data to good, bad and gray variants according to the statistics
    split_parser = subparsers.add_parser('split', help='variants splitting help')
    split_parser_required = split_parser.add_argument_group('required arguments')
    split_parser_required.add_argument('-i', '--input', dest='raw_stat_file', required=True, help='a raw variant dataset for splitting')
    split_parser.add_argument('-o', '--output', dest='output_handle', default=None, required=False,
                              help='the latter part of output filename [optional]')
    split_parser.add_argument('-m', '--model', dest='model', required=False, default='A',
                              help='splitting method for a specific random forest model, A or B [optional]. default: A')

    # arguments for variant classification
    classify_parser = subparsers.add_parser('classify', help='variants classification help')
    classify_parser_required = classify_parser.add_argument_group('required arguments')
    classify_parser_required.add_argument('-g', '--good', dest='good_var', required=True, help='good variants')
    classify_parser_required.add_argument('-b', '--bad', dest='bad_var', required=True, help='bad variants')
    classify_parser_required.add_argument('-y', '--gray', dest='gray_var', required=True, help='gray variants')
    classify_parser.add_argument('-m', '--model', dest='model', required=False, default='A',
                                 help='type of random forest model, A or B [optional]. default: A')
    classify_parser.add_argument('-o', '--output', dest='output_suffix', required=False, default='variants.tsv',
                                 help='the latter part of output filename [optional]. default: variants.tsv')

    # arguments for set outlier_gq and outlier_dp
    setOutlier_parser = subparsers.add_parser('setOutlier', help='set outlier help')
    setOutlier_parser_required = setOutlier_parser.add_argument_group('required arguments')
    setOutlier_parser_required.add_argument('-i', '--input', required=True, dest='file_list',
                                            help='input file(s). please separate the files with comma if there'
                                                 'are multiple files')

    args = vars(parser.parse_args())
    command = args.pop('command', None)
    return command, args

def main_setOutlier(**kwargs):
    file_list = kwargs['file_list'].split(',')
    setOutlier(file_list)

def main_stat(**kwargs):
    target_file = kwargs['target_file']
    stat_file = kwargs['stat_file']
    gc_file = kwargs['gc_file']
    ped_file = kwargs['ped_file']
    discord_geno_file = kwargs['discord_geno_file']
    hwe_file = kwargs['hwe_file']
    gender_file = kwargs['gender_info']

    global dp_threshold, gq_threshold
    dp_threshold = kwargs['dp']
    gq_threshold = kwargs['gq']

    # discordant genotype file processing
    discord_geno_dict = getDiscordInfo(discord_geno_file)
    # vcf file processing and data recording
    vcfProcessing(target_file, stat_file, gc_file, ped_file, discord_geno_dict, hwe_file, gender_file)

def main_split(**kwargs):
    input_file = kwargs['raw_stat_file']
    output_file = kwargs['output_handle']
    executeSplit(input_file, output_file)

def main_classify(**kwargs):
    good_var = kwargs['good_var']
    bad_var = kwargs['bad_var']
    gray_var = kwargs['gray_var']
    model = kwargs['model']
    output_suffix = kwargs['output_suffix']
    executeClassification(good_var, bad_var, gray_var, model, output_suffix)


def main():
    command_functions = {'stat': main_stat, 'split': main_split, 'classify': main_classify,
                         'setOutlier': main_setOutlier}
    command, args = parse_args()
    command_functions[command](**args)

if __name__ == '__main__':
    main()