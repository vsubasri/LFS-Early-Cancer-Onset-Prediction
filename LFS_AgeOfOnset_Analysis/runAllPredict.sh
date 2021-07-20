#!/bin/bash

module load R/3.6.1
seeds=(1342 2223 2345 5556 5678 6543 6667 7007 7778 8765 9009)

for seed in ${seeds[@]}
do

	outdir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/rds/Seed${seed}_TrainTestSplit_S25/
	models=$(find ${outdir}rds -maxdepth 1 -name "NoobCorrected*TrainingSet.rds")
	for model in ${models[@]}
	do 

		id=$(basename $model TrainingSet.rds)
		echo "[ Seed ]: $seed [ Running model ]: $id"
		covars="--sex --canceratdraw --systreat"
		
			echo "${id}"
			qsub   -v id=$id,vars="$covars --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --lfs --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --lfs --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --lfs --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --lfs --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --lfs --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --lfs --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
			qsub   -v id=$id,vars="$covars --lfs --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --cancer --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --cancer --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --cancer --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --cancer --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --cancer --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --cancer --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --cancer --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --lfs --cancer --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --lfs --cancer --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --lfs --cancer --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --lfs --cancer --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	        qsub   -v id=$id,vars="$covars --lfs --cancer --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --lfs --cancer --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --lfs --cancer --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --lfs --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --lfs --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --cancer --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --cancer --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --cancer --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --cancer --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --cancer --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --cancer --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --cancer --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --cancer --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --cancer --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --cancer --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --cancer --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --cancer --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --cancer --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
            qsub   -v id=$id,vars="$covars --scale --lfs --cancer --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
	done
done

