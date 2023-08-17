#!/bin/bash
#SBATCH --mem=118G
#SBATCH --time 48:00:00
#SBATCH --job-name bumphunter

module load R/3.4.4

Rscript /hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset/bumphunter.R


