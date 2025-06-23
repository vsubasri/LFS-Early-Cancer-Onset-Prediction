#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time 48:00:00
#SBATCH --job-name combine_seeds
#SBATCH -e log.combine.err
#SBATCH -o log.combine.out

module load R/3.6.1

seeds=(1 2 3 4 5 )
nsplit=40
scripts_dir="$(dirname "$0")"

for seed in ${seeds[@]}
do

	id=Seed${seed}_TrainTestSplit_S${nsplit}

	maindir="$(dirname "$0")/../data"
	dir=$maindir/$id/
	Rscript ${scripts_dir}/combine_results.R --id $id --dir $dir

done

Rscript ${scripts_dir}/get_average_results.R --dir $maindir

