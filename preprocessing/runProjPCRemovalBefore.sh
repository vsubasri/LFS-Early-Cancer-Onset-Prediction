#!/bin/bash
#SBATCH --mem=48G
#SBATCH --time 48:00:00
#SBATCH --job-name PCAreconstructionBefore
#SBATCH -e log.PCAreconstructionBefore.err
#SBATCH -o log.PCAreconstructionBefore.out

module load R/3.6.1

cd ${outdir} 

scripts_dir="$(dirname "$0")"

Rscript ${scripts_dir}/PCAreconstruction_Project.R \
        --id Noob_${methvalue} \
	--seed $seed \
	--outdir ${outdir} \
	--nsplit $nsplit

