#!/bin/bash
#SBATCH --mem=48G
#SBATCH --time 24:00:00
#SBATCH --job-name single_model
#SBATCH -e log.model.err
#SBATCH -o log.model.out

module load R/3.6.1

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

echo "Rscript ${scripts_dir}/runSingleModel.R --id $id $vars --outdir $outdir"

Rscript ${scripts_dir}/runSingleModel.R --id $id $vars --outdir $outdir

