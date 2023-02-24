#!/bin/bash

#SBATCH -t 200:00:00
#SBATCH -N 1 -c 40
#SBATCH --mem=50g
#SBATCH -o /hpf/largeprojects/davidm/blaverty/ageofonset/logs
#SBATCH -e /hpf/largeprojects/davidm/blaverty/ageofonset/logs


module load R/4.0.0

Rscript rgset.R
