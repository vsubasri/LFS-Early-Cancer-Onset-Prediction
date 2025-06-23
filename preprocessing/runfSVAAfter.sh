#!/bin/bash
#SBATCH --mem=48G
#SBATCH --time 48:00:00
#SBATCH --job-name fSVAAfter
#SBATCH -e log.fSVAAfter.err
#SBATCH -o log.fSVAAfter.out

module load R/3.6.1

cd ${outdir}

scripts_dir="$(dirname "$0")"

Rscript ${scripts_dir}/fSVA_Project.R \
	--id NoobCorrected_covs_${methvalue} \
	--outcome ${outcome} \
	--seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

Rscript ${scripts_dir}/fSVA_Project.R \
        --id NoobCorrected_${methvalue} \
	--outcome ${outcome} \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

