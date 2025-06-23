#!/bin/bash

#SBATCH -t 200:00:00
#SBATCH -N 1 -c 40
#SBATCH --mem=50g
#SBATCH -o logs/
#SBATCH -e logs/


module load R/4.0.0

Rscript rgset.R
