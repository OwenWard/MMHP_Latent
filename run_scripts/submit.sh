#!/bin/sh
#
#SBATCH -A stats
#SBATCH -J dsnl
#SBATCH -c 4
#SBATCH -t 320:00
#SBATCH --mem-per-cpu 8gb
#SBATCH -a 2-10
#SBATCH --mail-user=ogw2103@columbia.edu
#SBATCH --mail-type=ALL
##SBATCH -o cmmhp-%A.%a.out
##SBATCH --error=cmmhp-%A.%a.err 

module load R/4.0.1
echo "Launching R"
date

##R CMD BATCH --no-save --vanilla overall_predict.R routput_ov_pred_long$SLURM_ARRAY_TASK_ID
R CMD BATCH --no-save --vanilla active_dsnl.R routput_active_dsnl_$SLURM_ARRAY_TASK_ID
echo "Completed"
echo "Completed"
date

#end of script
