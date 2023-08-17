#!/bin/bash
#SBATCH --mem=72G
#SBATCH --time 48:00:00
#SBATCH --job-name BatchCorrectionBefore
#SBATCH -e log.BatchCorrectionBefore.err
#SBATCH -o log.BatchCorrectionBefore.out

module load R/3.6.1

cd ${outdir}

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

Rscript ${scripts_dir}/map450_850.R \
        --datafile ${outdir}rds/Noob_${methvalue}.rds \
        --value covs_${methvalue} \
        --covar T \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

Rscript ${scripts_dir}/map450_850.R \
        --datafile ${outdir}rds/Noob_${methvalue}.rds \
        --value ${methvalue} \
        --covar F \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

