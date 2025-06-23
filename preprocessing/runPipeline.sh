#!/bin/bash

methvalue=beta
seeds=(1 2 3 4 5)
##6 7 8 9 10)
nsplit=50
data_dir="$(dirname "$0")/../data"

for seed in ${seeds[@]}
do

	outdir=${data_dir}/LFS_ageofonset/LFS_with_COG/Seed${seed}_TrainTestSplit_S${nsplit}/
#	outdir=${data_dir}/LFS_ageofonset/new_clinical_COG_S${nsplit}/Seed${seed}_TrainTestSplit_S${nsplit}/

	mkdir ${outdir}
	mkdir ${outdir}rds
	mkdir ${outdir}Plots
	mkdir ${outdir}Output
	
	cp ${data_dir}/rds/Noob_beta.rds ${outdir}rds

	dep_bc=($(sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit runBatchCorrectionBefore.sh))
	dep_baPred=($(sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit runbaPredBefore.sh))
	dep_PC=($(sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit runProjPCRemovalBefore.sh))
	dep_fSVA_age=($(sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,outcome=ageofonset runfSVABefore.sh))
	dep_fSVA_tissue=($(sbatch --export methvalue=$methvalue,seed=$seed,nsplit=$nsplit,outdir=$outdir,outcome=tissue_type runfSVABefore.sh))
	dep_fSVA_cancer=($(sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,outcome=cancer_diagnosis runfSVABefore.sh))

	sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,id=ProjPC2Adj --dependency afterok:${dep_PC[3]} runBatchCorrectionAfter.sh
	sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,id=baPredComBat --dependency afterok:${dep_baPred[3]} runBatchCorrectionAfter.sh
	sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,id=fSVAageofonset --dependency afterok:${dep_fSVA_age[3]} runBatchCorrectionAfter.sh
	sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,id=fSVAtissue_type --dependency afterok:${dep_fSVA_tissue[3]} runBatchCorrectionAfter.sh
	sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,id=fSVAcancer_diagnosis --dependency afterok:${dep_fSVA_cancer[3]} runBatchCorrectionAfter.sh
	sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit --dependency afterok:${dep_bc[3]} runbaPredAfter.sh
	sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit --dependency afterok:${dep_bc[3]} runProjPCRemovalAfter.sh
	sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,outcome=ageofonset --dependency afterok:${dep_bc[3]} runfSVAAfter.sh
        sbatch --export methvalue=$methvalue,seed=$seed,outdir=$outdir,outcome=tissue_type --dependency afterok:${dep_bc[3]} runfSVAAfter.sh
        sbatch --export  methvalue=$methvalue,seed=$seed,outdir=$outdir,outcome=cancer_diagnosis --dependency afterok:${dep_bc[3]} runfSVAAfter.sh

done

