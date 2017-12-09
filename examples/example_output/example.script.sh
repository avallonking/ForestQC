#!/bin/bash

ForestQC set_outlier -i ../test.vcf.gz
ForestQC stat -i ../test.vcf.gz -o example.result.tsv -c ../gc_content_hg19.tsv -af ../user_features.tsv -p ../PedStructSeqID.txt --dp 20.0 --gq 51.0
ForestQC split -i example.result.tsv -af ../user_features.tsv -t ../user_thresholds.tsv 
ForestQC classify -g good.example.result.tsv -b bad.example.result.tsv -y gray.example.result.tsv -af ../user_features.tsv
