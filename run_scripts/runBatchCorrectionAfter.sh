#!/bin/bash
#SBATCH --mem=48G
#SBATCH --time 48:00:00
#SBATCH --job-name BatchCorrectionAfter
#SBATCH -e log.BatchCorrectionAfter.err
#SBATCH -o log.BatchCorrectionAfter.out

module load R/3.6.1

cd ${outdir}

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

Rscript ${scripts_dir}/map450_850.R \
        --datafile ${outdir}rds/Noob_${methvalue}_${id}.rds \
        --value after_covs_${methvalue}_${id} \
        --covar T \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

Rscript ${scripts_dir}/map450_850.R \
        --datafile ${outdir}rds/Noob_${methvalue}_${id}.rds \
        --value after_${methvalue}_${id} \
        --covar F \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

