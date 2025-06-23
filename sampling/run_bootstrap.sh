#!/bin/bash

id=NoobCorrected_after_covs_beta_ProjPC2Adj_
model=xgboost
outdir="$(dirname "$0")/../data/Seed6543_TrainTestSplit_S25/"
x=0
probes=220
genes=147
covars="--scale"
while [ $x -le 100 ]
do
    echo $x
    x=$(($x + 1))
    scripts_dir="$(dirname "$0")"
dep=$(qsub -v id=$id,bootstrap=$x,model=$model,outdir=$outdir,covars=${covars},probes=$probes,genes=$genes "${scripts_dir}/sample_refit.sh")
qsub -W depend=afterok:$dep -v id=$id,bootstrap=$x,model=$model,outdir=$outdir,covars=${covars},probes=$probes,genes=$genes "${scripts_dir}/sample_refit_predict.sh"

done



