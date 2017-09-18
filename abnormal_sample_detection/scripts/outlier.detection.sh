#$ -N outlier
#$ -cwd
#$ -o /u/scratch2/k/k8688933/github_repo/classifier/tests/sample_level_qc/log/outlier.detection.out
#$ -j y
#$ -m n
#$ -l h_data=10G,h_rt=5:00:00,highp
#$ -q eeskin_pod_16.q,eeskin_pod_12.q

outdir=/u/scratch2/k/k8688933/github_repo/classifier/tests/sample_level_qc
dp_file_list=/u/scratch2/k/k8688933/github_repo/classifier/tests/sample_level_qc/dp.file.list
abhet_file_list=/u/scratch2/k/k8688933/github_repo/classifier/tests/sample_level_qc/abhet.file.list
find_outlier_abhet_script=/u/scratch2/k/k8688933/github_repo/classifier/abnormal_sample_detection/scripts/find_outlier_abhet.R
find_outlier_depth_script=/u/scratch2/k/k8688933/github_repo/classifier/abnormal_sample_detection/scripts/find_outlier_dp.R

. /u/local/Modules/default/init/modules.sh
module load R/3.3.0

output=$1
picfile=plot.dp.$output\.png
output_dp=dp.$output
output_abhet=abhet.$output

ls $outdir/*.ab.csv > $abhet_file_list
ls $outdir/*.dp.csv > $dp_file_list

Rscript $find_outlier_depth_script $dp_file_list $output_dp $picfile
Rscript $find_outlier_abhet_script $abhet_file_list $output_abhet
