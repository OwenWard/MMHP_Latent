#!/bin/sh
#
#SBATCH -A stats
#SBATCH -J cmmhp
#SBATCH -c 4
#SBATCH -t 720:00
#SBATCH --mem-per-cpu 8gb
#SBATCH -a 1-10
#SBATCH --mail-user=ogw2103@columbia.edu
#SBATCH --mail-type=ALL
##SBATCH -o cmmhp-%A.%a.out
##SBATCH --error=cmmhp-%A.%a.err 

module load R
echo "Launching R"
date

R CMD BATCH --no-save --vanilla rank_stability_exp.R routput_rank_sim_$SLURM_ARRAY_TASK_ID
echo "Completed"
date

#end of script
