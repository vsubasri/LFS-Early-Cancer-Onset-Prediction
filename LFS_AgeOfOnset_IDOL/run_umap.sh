#!/bin/bash
#SBATCH --job-name="umap"
#SBATCH -t 10:00:00
#SBATCH -N 1 -c 30
#SBATCH --mem 50g
#SBATCH -o /hpf/largeprojects/davidm/blaverty/ageofonset/logs/umap_o
#SBATCH -e /hpf/largeprojects/davidm/blaverty/ageofonset/logs/umap_e


module load R/4.0.0

Rscript umap.R 
