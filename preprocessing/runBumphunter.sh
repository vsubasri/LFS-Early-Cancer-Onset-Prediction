#!/bin/bash
#SBATCH --mem=118G
#SBATCH --time 48:00:00
#SBATCH --job-name bumphunter

module load R/3.4.4

Rscript "$(dirname "$0")/../model/bumphunter.R"


