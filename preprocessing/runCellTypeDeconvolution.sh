#!/bin/bash
#SBATCH --mem=72G
#SBATCH --time 72:00:00
#SBATCH --job-name CellTypeDeconvolution

module purge; module load R/3.4.4

scripts_dir="$(dirname "$0")"

Rscript ${scripts_dir}/celltype_deconvolution.R


