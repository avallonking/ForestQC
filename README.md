# Classifiers

This classifier uses random forest model to identify good or bad variants from grey variants.
User can read the instruction in text format in doc/README

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

First, we need to calculate ths statistics from vcf file. **Note:**
 - The vcf file should have the information of each individual. 
 - We don't have to merge all vcf files together.
 - If no discordant genotype file provided, the number of discordant genotype of all variants will be NA
 - If no pedigree file provided, the mendel errors of all variants will be NA

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/stat.py -i [input_vcf] -o [output_filename] -p [ped_file(optional)] -d [discordant_genotype_file(optional)] -w [hwe_file(optional)]
```

You can check the usage of stat.py with

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/stat.py -h
```

Second, we need to divide the dataset into good, bad and grey variants. User can easily change the thresholds and output filename in the source. **Note that you don't have to merge the input file together**

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/data_preprocessing.py [input_file]
```

Third, we can train our random forest model and apply it on the classification. **Note that you need to merge all good variants into one file, all bad variants into one file and all grey variants into one file**

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/classification.py [good_variants] [bad_variants] [grey_variants] [output_filename_suffix]
```

### File format
 - Discordant genotype file (2 columns, tab separated)
  ```sh
  SNP_ID  Number_of_Discordant_Genotype
  chr1:259   1
  chr3:122   2
  ```

 - HWE p-value file (2 columns, tab separated):
  ```sh
  SNP_ID  HWE_p-value
  chr3:899   1.0
  cgr2:900   0.77
  ```
