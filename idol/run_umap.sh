#!/bin/bash
#SBATCH --job-name="umap"
#SBATCH -t 10:00:00
#SBATCH -N 1 -c 30
#SBATCH --mem 50g
#SBATCH -o logs/umap_o
#SBATCH -e logs/umap_e


module load R/4.0.0

Rscript umap.R 
