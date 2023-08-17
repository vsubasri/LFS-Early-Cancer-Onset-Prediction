#!/bin/bash
#SBATCH --mem=48G
#SBATCH --time 48:00:00
#SBATCH --job-name fSVABefore
#SBATCH -e log.fSVABefore.err
#SBATCH -o log.fSVABefore.out

module load R/3.6.1

cd ${outdir}

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

Rscript ${scripts_dir}/fSVA_Project.R \

	--id Noob_${methvalue} \
	--outcome ${outcome} \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

