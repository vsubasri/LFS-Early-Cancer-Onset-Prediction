#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.Bumphunter.err
#PBS -o log.Bumphunter.out
#PBS -l vmem=118g
#PBS -l walltime=48:00:00

module load R/3.4.4

cd /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runBumphunter.sh

