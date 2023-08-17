#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time 48:00:00
#SBATCH --job-name combine_seeds
#SBATCH -e log.combine.err
#SBATCH -o log.combine.out

module load R/3.6.1

seeds=(1 2 3 4 5 )
nsplit=40
scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

for seed in ${seeds[@]}
do

	id=Seed${seed}_TrainTestSplit_S${nsplit}
#	maindir=/hpf/largeprojects/davidm/vsubasri/methyl_data/LFS_ageofonset/LFS_with_COG
	maindir=/hpf/largeprojects/davidm/vsubasri/methyl_data/LFS_ageofonset/new_clinical_COG_S40
	dir=$maindir/$id/
	Rscript ${scripts_dir}/combine_results.R --id $id --dir $dir

done

Rscript ${scripts_dir}/get_average_results.R --dir $maindir

