#!/bin/bash
#PBS -S /bin/bash
#PBS -e Preprocessing.err
#PBS -o Preprocessing.out
#PBS -l vmem=72g
#PBS -l mem=72g
#PBS -l walltime=72:00:00

module load R/3.4.4

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/preprocess.R --value $methvalue

