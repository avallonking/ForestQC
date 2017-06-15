# Classifiers

This classifier uses random forest model to identify good or bad variants from grey variants.

### Requirements
 - Software: python > 3.3
 - Packages: scikit-learn, pandas, numpy

To install the packages
```sh
$ pip3 install -r $YOUR_PATH/classifier/requirements.txt
```

### Workflow
Before doing the classification, we should set the Outlier_GQ and Outlier_DP
```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/setOutlier.py [vcf_file1] [vcf_file2] [...]
```

Then go to vcf_stat.py, change DP_THRESHOLD and GQ_THRESHOLD to Outlier_DP and Outlier_GQ given by the previous step, respectively.

First, we need to calculate ths statistics from vcf file. Note the vcf file should have the information of each individual.

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/stat.py [vcf_file] [output_file] [ped_file] [discordant_genotype_file (optional)]
```

Second, we need to divide the dataset into good, bad and grey variants. User can easily change the thresholds and output filename in the source.
```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/data_preprocessing.py [input_file]
```
Third, we can train our random forest model and apply it on the classification
```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/classification.py [good_variants] [bad_variants] [grey_variants] [output_filename_suffix]
```
