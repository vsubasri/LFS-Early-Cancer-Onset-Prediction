#!/bin/bash
#SBATCH --mem=8G
#SBATCH --time 24:00:00
#SBATCH --job-name single_model
#SBATCH -e log.model.err
#SBATCH -o log.model.out

module load R/3.6.1

# Use relative path instead of absolute path
scripts_dir="$(dirname "$0")"
outdir=${outdir-"$(dirname "$0")/../data"}

echo "Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id $vars"
Rscript $scripts_dir/runSingleModel.R --outdir $outdir --id $id $vars

