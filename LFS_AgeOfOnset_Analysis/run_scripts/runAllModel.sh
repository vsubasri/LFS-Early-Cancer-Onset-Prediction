#!/bin/bash

nsplit=40
seeds=(1 2 3 4 5 6 7 8 9 10)

maindir=/hpf/largeprojects/davidm/vsubasri/methyl_data/LFS_ageofonset/new_clinical_COG_S40_auroc/
scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset
covars="--sex --tp53"
pred_covars="--fnw 2"

for seed in ${seeds[@]}
do
	outdir=${maindir}/Seed${seed}_TrainTestSplit_S${nsplit}/
	models=$(find ${outdir}rds -maxdepth 1 -name "NoobCorrected*ValidationSet.rds")
	for model in ${models[@]}
	do 

		id=$(basename $model ValidationSet.rds)
		echo "[ Seed ]: $seed [ Running model ]: $id"

		sbatch --export outdir=$outdir,covars="$covars",pred_covars="$pred_covars",id=$id ${scripts_dir}/runSeed.sh
	done
done

