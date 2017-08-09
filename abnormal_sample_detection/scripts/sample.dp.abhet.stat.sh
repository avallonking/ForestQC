#$ -N dp
#$ -cwd
#$ -o /u/scratch2/k/k8688933/github_repo/classifier/tests/sample_level_qc/log/dp.abhet.stat.out
#$ -j y
#$ -m n
#$ -t 1-305:1
# -l h_data=10G,h_rt=24:00:00,highp
#$ -l h_data=10G,h_rt=12:30:00,highp
#$ -q eeskin_pod_16.q,eeskin_pod_12.q

region=`expr $SGE_TASK_ID - 1`
window_size=10000
index_file=/u/scratch2/k/k8688933/github_repo/classifier/abnormal_sample_detection/support_files/indexfile
dir=/u/nobackup/eeskin2/jhsul/bipolar/gvcf_batch_final/snp
file_prefix=bp.churchill.hc.454indiv.snp.chr
outdir=/u/scratch2/k/k8688933/github_repo/classifier/tests/sample_level_qc/
depth_stat_script=/u/scratch2/k/k8688933/github_repo/classifier/abnormal_sample_detection/scripts/poolDP.py
abhet_stat_script=/u/scratch2/k/k8688933/github_repo/classifier/abnormal_sample_detection/scripts/calculate_sampleAB.py

i=-1
chr=0
temp=0
while IFS=' ', read xx yy; do
  i=`expr $i + 1`
  if [[ $i == $region ]]; then
    chr=$xx
    temp=$yy
  fi
done < $index_file

vcf_file=$dir/chr$chr/orig/$file_prefix.$chr.region.$region.vcf.gz
output_dp=$outdir/$chr.$region.dp.csv
output_ab=$outdir/$chr.$region.ab.csv

. /u/local/Modules/default/init/modules.sh
module load python/3.6.1

echo $vcf_file

python3 $depth_stat_script $vcf_file $output_dp $window_size
python3 $abhet_stat_script $vcf_file $output_ab
