#!/bin/bash
#SBATCH --mem=48G
#SBATCH --time 48:00:00
#SBATCH --job-name baPredAfter
#SBATCH -e log.baPredComBatAfter.err
#SBATCH -o log.baPredComBatAfter.out

module load R/3.6.1

cd ${outdir}

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

Rscript ${scripts_dir}/baPred_ComBat.R \
        --id NoobCorrected_covs_${methvalue} \
        --seed $seed \
	--outdir ${outdir} \
	--nsplit $nsplit

Rscript ${scripts_dir}/baPred_ComBat.R \
        --id NoobCorrected_${methvalue} \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

