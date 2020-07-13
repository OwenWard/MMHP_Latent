#!/bin/sh
#
#SBATCH -A tzsts
#SBATCH -J immhp
#SBATCH -c 8
#SBATCH -t 720:00
#SBATCH --mem-per-cpu 12gb
#SBATCH -a 1-10
#SBATCH --mail-user=ogw2103@columbia.edu
#SBATCH --mail-type=ALL

module load R
echo "Launching R"
date

R CMD BATCH --no-save --vanilla i_mmhp.R routput_immhp
echo "Completed"
date

#end of script
