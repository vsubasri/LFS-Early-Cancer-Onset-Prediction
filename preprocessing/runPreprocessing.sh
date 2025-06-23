#!/bin/bash
#SBATCH --mem=64G
#SBATCH --time 24:00:00
#SBATCH --job-name preprocessing
#SBATCH -e log.preprocessing.err
#SBATCH -o log.preprocessing.out

module load R/3.6.1

# Use relative path instead of absolute path
scripts_dir="$(dirname "$0")"

echo "[ Processing beta values ]"
Rscript ${scripts_dir}/preprocess_beta.R

echo "[ Processing M values ]"
Rscript ${scripts_dir}/preprocess_Mval.R

