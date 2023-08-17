#!/bin/bash
#SBATCH --mem=72G
#SBATCH --time 48:00:00
#SBATCH --job-name preprocessing
#SBATCH -e log.preprocessing.err
#SBATCH -o log.preprocessing.out

module load R/3.4.4

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

Rscript ${scripts_dir}/preprocess_beta.R

Rscript ${scripts_dir}/preprocess_Mval.R

