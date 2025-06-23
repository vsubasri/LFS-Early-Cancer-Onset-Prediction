#!/bin/bash
#SBATCH --mem=8G
#SBATCH --time 24:00:00
#SBATCH --job-name single_model
#SBATCH -e log.model.err
#SBATCH -o log.model.out

module load R/3.6.1

covars="--tp53"
script_dir="$(dirname "$0")"
outdir="$(dirname "$0")/../data"
id=NoobCorrected_beta_ProjPC2Adj_
region="3UTR"

echo "Rscript ${script_dir}/runSingleModel.R --outdir $outdir --id $id --lfs --aggregate $region $covars"
Rscript "${script_dir}/runSingleModel.R" --outdir $outdir --id $id --lfs --aggregate $region $covars

echo "Rscript ${script_dir}/predictSingleExtVal.R --outdir $outdir --id $id --lfs --aggregate $region $covars --fnw 1"
Rscript "${script_dir}/predictSingleExtVal.R" --outdir $outdir --id $id --lfs --aggregate $region $covars --fnw 1

