#!/bin/sh
#
#SBATCH -A tzsts
#SBATCH -J cmmhp
#SBATCH -c 8
#SBATCH -t 720:00
#SBATCH --mem-per-cpu 12gb
#SBATCH -a 1-10

module load R
echo "Launching R"
date

R CMD BATCH --no-save --vanilla c_mmhp.R routput_cmmhp
echo "Completed"
date

#end of script
