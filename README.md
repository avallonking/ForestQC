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
Before doing the classification, we should set the Outlier_GQ and Outlier_DP. It would print out Outlier_DP and Outlier_GQ on the screen, which would be used as inputs in the next step. *Note that all vcf files in this analysis should be included.*
```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/setOutlier.py [vcf_file1] [vcf_file2] [...]
```

Then go to vcf_stat.py, set Outlier_DP and Outlier_GQ given by the previous step, respectively.

First, we need to calculate ths statistics from vcf file. *It will output a file containing all statistics information for each variant.* **Note:**
 - The vcf file should have the information of each individual. 
 - We don't have to merge all vcf files together.
 - If no discordant genotype file provided, the number of discordant genotype of all variants will be NA
 - If no pedigree file provided, the mendel errors of all variants will be NA

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/stat.py -i [input_vcf] -o [output_filename] -g [gender_file(optional)] -p [ped_file(optional)] -d [discordant_genotype_file(optional)] -w [hwe_file(optional)] --gq [Outlier_GQ] --dp [Outlier_DP]
```

You can check the usage of stat.py with

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/stat.py -h
```

Second, we need to divide the dataset into good, bad and grey variants. User can easily change the thresholds and output filename in the source. *The output files would be three: good.xx, bad.xx and grey.xx.* **Note that you don't have to merge the input file together**

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/data_preprocessing.py [input_file] [output_filename_suffix (optional)]
```

Third, we can train our random forest model and apply it on the classification of grey variants. *The output files would be predicted_good_xx and predictred_bad_xx.* **Note that you need to merge all good variants into one file, all bad variants into one file and all grey variants into one file**

```sh
$ python3 $YOUR_PATH/classifier/random_forest_classifier/classification.py [good_variants] [bad_variants] [grey_variants] [output_filename_suffix]
```

###File format
 - Gender file(2 columns, tab separated. f for female, m for male):
 ```sh
 SampleID  Gender
 12312sd12  m
 4342sd423  f
 ```

 - Pedigree file(9 columns, tab separated)
 ```sh
 FamilyID IndividualID FatherID MotherID Sex PhenotypeID(1 if it is control) DBPID SampleID SeqID
 C1  2 3 4 m  1  xxx xxx xxx
 ```

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

 - Result file (The variant is a good variant if Good = 1, or a bad variant if Good = 0, grey variants do not have Good column before it is predicted)
  ```sh
  RSID CHR POS REF ALT MAF Mean_DP Mean_GQ SD_DP SD_GQ Outlier_DP Outlier_GQ Discordant_Geno Mendel_Error Missing_Rate HWE ABHet ABHom Good
  chr1:144  1 144 A T 0.03  54.00 54.00 23.00 13.24 0.43  0.23  1 3 0.01  1.0 0.45  0.99  1
  ```
