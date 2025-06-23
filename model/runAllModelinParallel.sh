#!/bin/bash

module load R/3.6.1
seeds=(1 2 3 4 5)
nsplit=30
script_dir="$(dirname "$0")"
covars="--canceratdraw --systreat --scale"
covars_label="scaled_canceratdraw_systreat"

for seed in ${seeds[@]}
do
	outdir="$(dirname "$0")/../data/Seed${seed}_TrainTestSplit_S$nsplit/"
	mkdir -p ${outdir}rds
	models=$(find ${outdir}rds -maxdepth 1 -name "NoobCorrected*TrainingSet.rds" 2>/dev/null || echo "")
	if [ -z "$models" ]; then
		echo "No models found for seed $seed, skipping"
		continue
	fi
	
	for model in ${models[@]}
	do 

		id=$(basename $model TrainingSet.rds)
		echo "[ Seed ]: $seed [ Running model ]: $id"
		
                if [ ! -f ${outdir}rds/${id}gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			echo "${id}gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds"
	        	dep=($(sbatch --export id=$id,vars="$covars --gene",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --gene",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi

                if [ ! -f ${outdir}rds/${id}Body_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=($(sbatch --export id=$id,vars="$covars --aggregate Body",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --aggregate Body",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi

                if [ ! -f ${outdir}rds/${id}TSS200_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=($(sbatch --export id=$id,vars="$covars --aggregate TSS200",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --aggregate TSS200",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}TSS1500_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=($(sbatch --export id=$id,vars="$covars --aggregate TSS1500",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --aggregate TSS1500",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}3UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=($(qsub --export id=$id,vars="$covars --aggregate 3UTR",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --aggregate 3UTR",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}5UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=($(sbatch --export id=$id,vars="$covars --aggregate 5UTR",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --aggregate 5UTR",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}1stExon_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=($(sbatch --export id=$id,vars="$covars --aggregate 1stExon",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --aggregate 1stExon",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi

		if [ ! -f ${outdir}rds/${id}lfs_TSS200_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=($(qsub --export id=$id,vars="$covars --lfs --aggregate TSS200",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --aggregate TSS200",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
                if [ ! -f ${outdir}rds/${id}lfs_TSS1500_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=($(qsub --export id=$id,vars="$covars --lfs --aggregate TSS1500",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --aggregate TSS1500",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
                if [ ! -f ${outdir}rds/${id}lfs_3UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=($(sbatch --export id=$id,vars="$covars --lfs --aggregate 3UTR",outdir=$outdir ${script_dir}/runSingleModel.sh))
                        sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --aggregate 3UTR",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
                if [ ! -f ${outdir}rds/${id}lfs_5UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=($(sbatch --export id=$id,vars="$covars --lfs --aggregate 5UTR",outdir=$outdir ${script_dir}/runSingleModel.sh))
                        sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --aggregate 5UTR",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}lfs_1stExon_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
		        dep=($(sbatch --export id=$id,vars="$covars --lfs --aggregate 1stExon",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --aggregate 1stExon",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}lfs_gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
		        dep=($(qsub --export id=$id,vars="$covars --lfs --gene",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --gene",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}lfs_Body_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
		        dep=($(sbatch --export id=$id,vars="$covars --lfs --aggregate Body",outdir=$outdir ${script_dir}/runSingleModel.sh))
			sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --aggregate Body",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi

		if [ ! -f ${outdir}rds/${id}nocancer_TSS200_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --cancer --aggregate TSS200",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --cancer --aggregate TSS200",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}nocancer_TSS1500_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --cancer --aggregate TSS1500",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --cancer --aggregate TSS1500",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}nocancer_1stExon_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --cancer --aggregate 1stExon",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --cancer --aggregate 1stExon",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi             
		if [ ! -f ${outdir}rds/${id}nocancer_gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --cancer --gene",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --cancer --gene",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}nocancer_Body_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
               	dep=($(sbatch --export id=$id,vars="$covars --cancer --aggregate Body",outdir=$outdir ${script_dir}/runSingleModel.sh))
               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --cancer --aggregate Body",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
                if [ ! -f ${outdir}rds/${id}nocancer_3UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=($(sbatch --export id=$id,vars="$covars --cancer --aggregate 3UTR",outdir=$outdir ${script_dir}/runSingleModel.sh))
                        sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --cancer --aggregate 3UTR",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
                if [ ! -f ${outdir}rds/${id}nocancer_5UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=($(sbatch --export id=$id,vars="$covars --cancer --aggregate 5UTR",outdir=$outdir ${script_dir}/runSingleModel.sh))
                        sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --cancer --aggregate 5UTR",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi

		if [ ! -f ${outdir}rds/${id}nocancer_lfs_TSS200_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --lfs --cancer --aggregate TSS200",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --cancer --aggregate TSS200",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}nocancer_lfs_TSS1500_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --lfs --cancer --aggregate TSS1500",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --cancer --aggregate TSS1500",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}nocancer_lfs_1stExon_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --lfs --cancer --aggregate 1stExon",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --cancer --aggregate 1stExon",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}nocancer_lfs_gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --lfs --cancer --gene",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --cancer --gene",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}nocancer_lfs_Body_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=($(sbatch --export id=$id,vars="$covars --lfs --cancer --aggregate Body",outdir=$outdir ${script_dir}/runSingleModel.sh))
	               	sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --cancer --aggregate Body",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
                if [ ! -f ${outdir}rds/${id}nocancer_lfs_3UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=($(sbatch --export id=$id,vars="$covars --lfs --cancer --aggregate 3UTR",outdir=$outdir ${script_dir}/runSingleModel.sh))
                        sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --cancer --aggregate 3UTR",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
                fi
                if [ ! -f ${outdir}rds/${id}nocancer_lfs_5UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=($(sbatch --export id=$id,vars="$covars --lfs --cancer --aggregate 5UTR",outdir=$outdir ${script_dir}/runSingleModel.sh))
                        sbatch --dependency afterok:${dep[3]} --export id=$id,vars="$covars --lfs --cancer --aggregate 5UTR",outdir=$outdir ${script_dir}/predictSingleExtVal.sh
		fi
	done
done

