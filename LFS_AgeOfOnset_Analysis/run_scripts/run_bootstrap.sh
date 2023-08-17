#!/bin/bash

id=NoobCorrected_after_covs_beta_ProjPC2Adj_
model=xgboost
outdir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/rds/Seed6543_TrainTestSplit_S25/
x=0
probes=220
genes=147
covars="--scale"
while [ $x -le 100 ]
do
    echo $x
    x=$(($x + 1))
    dep=$(qsub -v id=$id,bootstrap=$x,model=$model,outdir=$outdir,covars=${covars},probes=$probes,genes=$genes /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit.sh)
    qsub -W depend=afterok:$dep -v id=$id,bootstrap=$x,model=$model,outdir=$outdir,covars=${covars},probes=$probes,genes=$genes /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit_predict.sh

done



