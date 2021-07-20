#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.err
#PBS -o log.out
#PBS -l vmem=24g
#PBS -l mem=24g
#PBS -l walltime=48:00:00

module unload R ; module load R/3.6.1

cd /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/

seeds=(1234 1342 6543 7007 9009)
##40_seeds
#seeds=(1234 1342 2457 4567 6543 7007 8613 8642 9009 9876)
#25_seeds
#seeds=(1342 2223 2345 5556 5678 6543 6667 7007 7778 8765 9009)

for seed in ${seeds[@]}
do

	id=Seed${seed}_TrainTestSplit_S40
	dir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/rds/$id/
	Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/combineResults.R --id $id --dir $dir

done

