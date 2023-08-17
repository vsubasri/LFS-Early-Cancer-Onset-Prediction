#!/bin/bash
#SBATCH --mem=8G
#SBATCH --time 24:00:00
#SBATCH --job-name single_model
#SBATCH -e log.model.err
#SBATCH -o log.model.out

module load R/3.6.1

covars="--tp53"
scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset
outdir=/hpf/largeprojects/davidm/vsubasri/methyl_data/LFS_ageofonset/new_clinical_COG_S40/Seed4_TrainTestSplit_S40/
id=NoobCorrected_beta_ProjPC2Adj_
region="3UTR"

echo "Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --aggregate $region $covars"
Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id --lfs --aggregate $region $covars

echo "Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate $region $covars --fnw 1"
Rscript $scripts_dir/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate $region $covars --fnw 1

