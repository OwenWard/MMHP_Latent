#!/bin/sh
#
#SBATCH -A tzsts
#SBATCH -J cmmhp
#SBATCH -c 8
#SBATCH -t 1500:00
#SBATCH --mem-per-cpu 12gb
#SBATCH -a 1-10
#SBATCH --mail-user=ogw2103@columbia.edu
#SBATCH --mail-type=ALL

module load R
echo "Launching R"
date

R CMD BATCH --no-save --vanilla c_mmhp.R routput_cmmhp_dc
echo "Completed"
date

#end of script
