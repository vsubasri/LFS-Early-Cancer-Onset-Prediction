#!/bin/bash
#SBATCH --mem=28G
#SBATCH --time 24:00:00
#SBATCH --job-name predict
#SBATCH -e log.predict.err
#SBATCH -o log.predict.out

module load R/3.6.1

# Use the directory where this script is located
script_dir="$(dirname "$0")"

Rscript "${script_dir}/predictSingleExtVal.R" --id $id $vars --outdir $outdir


