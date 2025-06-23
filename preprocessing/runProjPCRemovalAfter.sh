#!/bin/bash
#SBATCH --mem=48G
#SBATCH --time 48:00:00
#SBATCH --job-name PCAreconstruction
#SBATCH -e log.PCAreconstructionAfter.err
#SBATCH -o log.PCAreconstructionAfter.out

module load R/3.6.1

cd ${outdir}

scripts_dir="$(dirname "$0")"

Rscript ${scripts_dir}/PCAreconstruction_Project.R \
        --id NoobCorrected_covs_${methvalue} \
        --seed $seed \
	--outdir $outdir \
	--nsplit $nsplit

Rscript ${scripts_dir}/PCAreconstruction_Project.R \
        --id NoobCorrected_${methvalue} \
        --seed $seed \
	--outdir $outdir \
	--nsplit $nsplit

