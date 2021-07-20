#!/bin/bash

methvalue=beta
seeds=(1342 2223 2345 5556 5678 6543 6667 7007 7778 8765 9009)
nsplit=25

for seed in ${seeds[@]}
do

	outdir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/rds/Seed${seed}_TrainTestSplit_S${nsplit}/

	mkdir $outdir ${outdir}rds ${outdir}Plots ${outdir}Output

	dep_preprocess=$(qsub -v methvalue=$methvalue runPreprocessing.sh)

	dep_bc=$(qsub -W depend=afterok:$dep_preprocess -v methvalue=$methvalue,seed=$seed,outdir=$outdir runBatchCorrectionBefore.sh)
	dep_baPred=$(qsub -W depend=afterok:$dep_preprocess -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit runbaPredBefore.sh)
	dep_PC=$(qsub -W depend=afterok:$dep_preprocess -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit runProjPCRemovalBefore.sh)
	dep_fSVA_age=$(qsub -W depend=afterok:$dep_preprocess -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,outcome=ageofonset runfSVABefore.sh)
	dep_fSVA_tissue=$(qsub -W depend=afterok:$dep_preprocess -v methvalue=$methvalue,seed=$seed,outdir=$outdir,outcome=tissue_type runfSVABefore.sh)
	dep_fSVA_cancer=$(qsub -W depend=afterok:$dep_preprocess -v methvalue=$methvalue,seed=$seed,outdir=$outdir,outcome=cancer_diagnosis runfSVABefore.sh)

	qsub -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,id=ProjPC2Adj -W depend=afterok:$dep_PC runBatchCorrectionAfter.sh
	qsub -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,id=baPredComBat -W depend=afterok:$dep_baPred runBatchCorrectionAfter.sh
	qsub -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,id=fSVAageofonset -W depend=afterok:$dep_fSVA_age runBatchCorrectionAfter.sh
	qsub -v methvalue=$methvalue,seed=$seed,outdir=$outdir,id=fSVAtissue_type -W depend=afterok:$dep_fSVA_tissue runBatchCorrectionAfter.sh
	qsub -v methvalue=$methvalue,seed=$seed,outdir=$outdir,id=fSVAcancer_diagnosis -W depend=afterok:$dep_fSVA_cancer runBatchCorrectionAfter.sh

	qsub -W depend=afterok:$dep_bc -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit runbaPredAfter.sh
	qsub -W depend=afterok:$dep_bc -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit runProjPCRemovalAfter.sh
	qsub -W depend=afterok:$dep_bc -v methvalue=$methvalue,seed=$seed,outdir=$outdir,nsplit=$nsplit,outcome=ageofonset runfSVAAfter.sh
    qsub -W depend=afterok:$dep_bc -v methvalue=$methvalue,seed=$seed,outdir=$outdir,outcome=tissue_type runfSVAAfter.sh
    qsub -W depend=afterok:$dep_bc -v methvalue=$methvalue,seed=$seed,outdir=$outdir,outcome=cancer_diagnosis runfSVAAfter.sh

done

