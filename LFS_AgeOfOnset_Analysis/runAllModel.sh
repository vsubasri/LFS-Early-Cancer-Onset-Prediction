#!/bin/bash


module load R/3.6.1

seeds=(1342 2223 2345 5556 5678 6543 6667 7007 7778 8765 9009)

lfs=true
sex=true
canceratdraw=true
systreat=true
scale=true

for seed in ${seeds[@]}
do

	outdir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/rds/Seed${seed}_TrainTestSplit_S25/
	models=$(find ${outdir}rds -maxdepth 1 -name "NoobCorrected*TrainingSet.rds")
	for model in ${models[@]}
	do 

		if [ "$lfs" = true ] ; then
    			covars="${covars} --lfs" 
    			covars_label="${covars}_lfs" 
		fi

		if [ "$sex" = true ] ; then
    			covars="${covars} --sex" 
    			covars_label="${covars}_sex" 
		fi

		if [ "$canceratdraw" = true ] ; then
    			covars="${covars} --canceratdraw" 
    			covars_label="${covars}_canceratdraw" 
		fi

		if [ "$systreat" = true ] ; then
    			covars="${covars} --systreat" 
    			covars_label="${covars}_systreat" 
		fi

		if [ "$scale" = true ] ; then
    			covars="${covars} --scale" 
    			covars_label="${covars}_scale"
		fi

		id=$(basename $model TrainingSet.rds)
		echo "[ Seed ]: $seed [ Running model ]: $id"
		covars="--sex --canceratdraw --systreat --scale"
		covars_label="scaled_sex_canceratdraw_systreat"
		
                if [ ! -f ${outdir}rds/${id}gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			echo "${id}gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds"
	        	dep=$(qsub -v id=$id,vars="$covars --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi

                if [ ! -f ${outdir}rds/${id}Body_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=$(qsub -v id=$id,vars="$covars --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi

                if [ ! -f ${outdir}rds/${id}TSS200_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=$(qsub -v id=$id,vars="$covars --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}TSS1500_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=$(qsub -v id=$id,vars="$covars --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}3UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=$(qsub -v id=$id,vars="$covars --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}5UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=$(qsub -v id=$id,vars="$covars --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}1stExon_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=$(qsub -v id=$id,vars="$covars --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi

		if [ ! -f ${outdir}rds/${id}lfs_TSS200_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=$(qsub -v id=$id,vars="$covars --lfs --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
                if [ ! -f ${outdir}rds/${id}lfs_TSS1500_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
			dep=$(qsub -v id=$id,vars="$covars --lfs --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
                if [ ! -f ${outdir}rds/${id}lfs_3UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=$(qsub -v id=$id,vars="$covars --lfs --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
                        qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
                if [ ! -f ${outdir}rds/${id}lfs_5UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=$(qsub -v id=$id,vars="$covars --lfs --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
                        qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}lfs_1stExon_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
		        dep=$(qsub -v id=$id,vars="$covars --lfs --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}lfs_gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
		        dep=$(qsub -v id=$id,vars="$covars --lfs --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}lfs_Body_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
		        dep=$(qsub -v id=$id,vars="$covars --lfs --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
			qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi

		if [ ! -f ${outdir}rds/${id}nocancer_TSS200_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --cancer --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --cancer --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}nocancer_TSS1500_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --cancer --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --cancer --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}nocancer_1stExon_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --cancer --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --cancer --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi             
		if [ ! -f ${outdir}rds/${id}nocancer_gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --cancer --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --cancer --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}nocancer_Body_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
               		dep=$(qsub -v id=$id,vars="$covars --cancer --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
               		qsub -W depend=afterok:$dep -v id=$id,vars="$covars --cancer --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
                if [ ! -f ${outdir}rds/${id}nocancer_3UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=$(qsub -v id=$id,vars="$covars --cancer --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
                        qsub -W depend=afterok:$dep -v id=$id,vars="$covars --cancer --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
                if [ ! -f ${outdir}rds/${id}nocancer_5UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=$(qsub -v id=$id,vars="$covars --cancer --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
                        qsub -W depend=afterok:$dep -v id=$id,vars="$covars --cancer --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi

		if [ ! -f ${outdir}rds/${id}nocancer_lfs_TSS200_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --lfs --cancer --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --cancer --aggregate TSS200",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
		if [ ! -f ${outdir}rds/${id}nocancer_lfs_TSS1500_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --lfs --cancer --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --cancer --aggregate TSS1500",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}nocancer_lfs_1stExon_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --lfs --cancer --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --cancer --aggregate 1stExon",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}nocancer_lfs_gene_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --lfs --cancer --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --cancer --gene",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
		if [ ! -f ${outdir}rds/${id}nocancer_lfs_Body_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
	               	dep=$(qsub -v id=$id,vars="$covars --lfs --cancer --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
	               	qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --cancer --aggregate Body",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi
                if [ ! -f ${outdir}rds/${id}nocancer_lfs_3UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=$(qsub -v id=$id,vars="$covars --lfs --cancer --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
                        qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --cancer --aggregate 3UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
                fi
                if [ ! -f ${outdir}rds/${id}nocancer_lfs_5UTR_${covars_label}_xgboost_ageofonset_ROCInfoVal.rds ] ; then
                        dep=$(qsub -v id=$id,vars="$covars --lfs --cancer --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.sh)
                        qsub -W depend=afterok:$dep -v id=$id,vars="$covars --lfs --cancer --aggregate 5UTR",outdir=$outdir /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.sh
		fi

	done
done

