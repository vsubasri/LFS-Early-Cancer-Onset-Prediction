#!/bin/bash
#SBATCH --mem=72G
#SBATCH --time 72:00:00
#SBATCH --job-name CellTypeDeconvolution

module purge; module load R/3.4.4

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

Rscript ${scripts_dir}/celltype_deconvolution.R


