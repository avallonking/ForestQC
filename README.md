# ForestQC

This classifier uses random forest model to identify good or bad variants from grey variants.

### Installation
To install the software. **This soft is compatible with OSX, 64-bit Linux and 64-bit Windows systems.**
```sh
$ conda install forestqc -c avallonking
```

After installation, you can check the usage 
```sh
$ ForestQC -h
# and check the usage for each command, for example, stat
$ ForestQC stat -h
```

### Workflow
```example/example_output/example.script.sh``` provides step-by-step tutorial about how to run this software.

Before doing the classification, we should set the Outlier_GQ and Outlier_DP. It would print out Outlier_DP and Outlier_GQ on the screen, which would be used as inputs in the next step. *Note that all vcf files in this analysis should be included.*
```sh
$ ForestQC set_outlier -i [vcf_file1] [vcf_file2] [...]
```

If you use other reference genome other than hg19, please use our script to generate a GC content table for your reference genome. Otherwise, we have prepared hg19 GC content table in ```examples/gc_content_hg19.tsv``` and you can directly use it.
```sh
# generate GC content table for your reference genome
$ ForestQC compute_gc -i [ref_genome.fasta] -o [output_file]
```

First, we need to calculate ths statistics from vcf file. *It will output a file containing all statistics information for each variant.* **Note:**
 - The vcf file should have the information of each individual. 
 - We don't have to merge all vcf files together.
 - If no discordant genotype file provided, the number of discordant genotype of all variants will be NA
 - If no pedigree file provided, the mendel errors of all variants will be NA

```sh
$ ForestQC stat -i [input_vcf] -o [output_filename] -c [gc_content_file] -g [gender_file(optional)] -p [ped_file (optional)] -d [discordant_genotype_file (optional)] -w [hwe_file(optional)] --gq [Outlier_GQ] --dp [Outlier_DP] -af [user_defined_features_file (optional)]
```

Second, we need to divide the dataset into good, bad and grey variants. User can easily change the thresholds and output filename in the source. *The output files would be three: good.xx, bad.xx and grey.xx.* **Note that you don't have to merge the input file together**

```sh
$ ForestQC split -i [input_file] -o [output_filename_suffix (optional)] -t [user_defined_threshold_file (optional)] -af [user_defined_features_names (if user-defined features added in last step, this is required)]
```

Third, we can train our random forest model and apply it on the classification of grey variants. *The output files would be predicted_good_xx and predictred_bad_xx.* **Note that you need to merge all good variants into one file, all bad variants into one file and all grey variants into one file**

```sh
$ ForestQC classify -g [good_variants] -b [bad_variants] -y [grey_variants] -o [output_filename_suffix (optional)] -af [user_defined_features_names (optional)]
```

### File format
#### Output file
 - Statistics file (tab-separated, assuming that we included user-defined features)
```sh
RSID CHR POS REF ALT MAF Mean_DP Mean_GQ SD_DP SD_GQ Outlier_DP Outlier_GQ Discordant_Geno Mendel_Error Missing_Rate HWE ABHet ABHom GC user-defined_feature1 user-defined_feature2
chr1:144  1 144 A T 0.03  54.00 54.00 23.00 13.24 0.43 0.23 1 3 0.01 1.0 0.45 0.99 0.435 2.0 4.3
chr1:144  1 144 A T 0.03  54.00 54.00 23.00 13.24 0.43 0.23 1 3 0.01 1.0 0.45 0.99 0.435 NA 4.3
```

 - Final result file (The variant is a good variant if Good = 1, or a bad variant if Good = 0, grey variants do not have Good column before it is predicted)
```sh
RSID CHR POS REF ALT MAF Mean_DP Mean_GQ SD_DP SD_GQ Outlier_DP Outlier_GQ Discordant_Geno Mendel_Error Missing_Rate HWE ABHet ABHom GC user-defined_feature1 user-defined_feature2 Good
chr1:144  1 144 A T 0.03  54.00 54.00 23.00 13.24 0.43 0.23 1 3 0.01 1.0 0.45 0.99 0.435 2.0 4.3 1
```
 
#### Input file
 - Gender file(2 columns, tab-separated. f for female, m for male):
```sh
SampleID  Gender
12312sd12  m
4342sd423  f
7hi987989  0     # Any samples whose gender is neither m nor f will be ignored
```

 - Pedigree file(9 columns, tab-separated)
```sh
FamilyID IndividualID FatherID MotherID Sex PhenotypeID(1 if it is control) DBPID SampleID SeqID
C1  2 3 4 m  1  xxx xxx xxx
```

 - Discordant genotype file (2 columns, tab-separated)
```sh
SNP_ID  Number_of_Discordant_Genotype
chr1:259   1
chr3:122   2
```

 - HWE p-value file (2 columns, tab-separated):
```sh
SNP_ID  HWE_p-value
chr3:899   1.0
cgr2:900   0.77
```

- User-defined features (tab-separated, data are calculated by users in advance. In practice, missing values will be imputed with the median of the sample)
```sh
RSID feature1 feature2 feature3
chr1:32 2.3 4.5 2.3
chr2:34 4.3 NA 5.6
```

- User-defined filter thresholds (tab-separated, 4 columns, "rare" means filters for rare variants, "common" for common variants, "all" for both)
```sh
good all filter threshold
bad all MAF threshold
bad rare Mendel_Error threshold
bad common Missing_Rate threshold
outlier rare Missing_Rate threshold
outlier common Missing_Rate threshold
```

- GC content file (tab-separated, 3 columns, including chrosome, starting positions and GC content)
```sh
chrM 1 0.3
chrM 1001 0.46
```
  
<!--
# Sample-level QC

This approach can detect samples with abnormal mean and standard deviation of sequencing depth, basd on the linear relationship between them. It can also find sample with extreme value of ABHet (0.4 <= ABHet <= 0.6 is consider as normal in default).

### Requirements
 - Software: python > 3.3, R >= 3.3.0
 - Packages(python): pandas, numpy
 - Packages(R): car, data.table

To install python packages:
```sh
$ pip3 install -r $YOUR_PATH/classifier/requirements.txt
```

To install R packages:
```sh
$ R
>>> install.packages(c("car", "data.table"))
```

### Usage
Sample usage is shown in ```./abnormal_sample_detection/scripts/sample.dp.abhet.stat.sh``` and ```./abnormal_sample_detection/scripts/outlier.detection.sh```. For simple usage, user can directly use these shell script for analysis on Hoffman2, but make sure to change some specific variables for your purpose. *Change the variables in these shell scripts:*
  - ```./abnormal_sample_detection/scripts/sample.dp.abhet.stat.sh```:
    - PATH for Hoffman2 log file (line 3)
    - Job array (line 6)
    - ```$indexfile``` (line 13)
    - ```$dir``` for vcf files (line 14)
    - ```$file_prefix``` for the prefix shared by all vcf files (line 15)
    - ```$outdir``` for the directory to store output files (line 16)
    - ```$depth_stat_script``` for the PATH of '''poolDP.py''' (line 17)
    - ```$abhet_stat_script``` for the PATH of '''calculate_sampleAB.py''' (line 18)
  - ```./abnormal_sample_detection/scripts/outlier.detection.sh```:
    - PATH for Hoffman2 log file (line 3)
    - ```$outdir``` for output directory
    - ```$dp_file_list``` and ```$abhet_file_list``` for the file list of all intermediate output files (line 10 and 11)
    - ```$find_outlier_abhet_script``` and ```$find_outlier_depth_script``` for the PATH of '''sh find_outlier_abhet.R''' and '''sh find_outlier_dp.R'''

For normal usage, please follow the 2 steps below:
  - Step1: Pool variants depth together into sliding windows of fixed size. And extract allele information from vcf files. **It is better to ignore chromosome X**.
  ```sh
  $ python3 ./abnormal_sample_detection/scripts/poolDP.py [vcf_file] [output_file] [window_size]
  $ python3 ./abnormal_sample_detection/scripts/calculate_sampleAB.py [vcf_file] [output_file]
  ```
  - Step2: Find outlier samples from mean and sd of depth, and ABHet.
  ```sh
  $ Rscript ./abnormal_sample_detection/scripts/find_outlier_dp.R [dp_file_list] [output_file] [png_file]
  $ Rscript ./abnormal_sample_detection/scripts/find_outlier_abhet.R [abhet_file_list] [output_file]
  ```

### Output
It will have temporary output and final output:
 - Temporary output:
   CSV files with depth information and CSV files with alllele counts.
 - Final ouput:
   TSV files with Mean, SD of depth and ABHet for each sample (need to concatenate output from 2 R programs with command paste). And a png file with a linear regression plot of Mean and SD, where outliers are labeled as red points.

### File format
  - file list (files are separated by \n)
  ```sh
  file1
  file2
  file3
  ```

  - depth file (comma separated, region_size is the number of variants in that region):
  ```sh
  ,sample1,sample2,sample3,region_size
  1,34.5,24.1,23.3,2423
  2,13,23,14.4,333
  ```

  - abhet file (comma separated):
  ```sh
  ,sample1,sample2,sample3
  ALT,325,2235,235
  REF,2342,1414,355
  ```

  - final output (tab separated):
  ```sh
           Mean  SD  ABHet
  sample1  24     8   0.4
  sample2  34     2   0.5
  ```
-->
