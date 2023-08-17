#!/bin/bash
#SBATCH --mem=28G
#SBATCH --time 24:00:00
#SBATCH --job-name predict
#SBATCH -e log.predict.err
#SBATCH -o log.predict.out

module load R/3.6.1

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

Rscript ${scripts_dir}/predictSingleExtVal.R --id $id $vars --outdir $outdir


