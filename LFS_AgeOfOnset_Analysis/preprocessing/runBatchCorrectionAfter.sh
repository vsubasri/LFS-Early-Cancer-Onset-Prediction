#!/bin/bash
#PBS -S /bin/bash
#PBS -e log/log.BatchCorrection.err
#PBS -o log/log.BatchCorrection.out
#PBS -l vmem=48g
#PBS -l mem=48g
#PBS -l walltime=48:00:00

export R_LD_LIBRARY_PATH=/hpf/tools/centos6/libgsl/2.1/lib
export R_LIBS_USER=/hpf/largeprojects/adam/valli/R/3.4.4/
export R_LIBS=/hpf/largeprojects/adam/valli/R/3.4.4/

module load R/3.4.4

cd /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/map450_850.R \
        --datafile ${outdir}rds/Noob_${methvalue}_${id}.rds \
        --id after_covs_${methvalue}_${id} \
        --covar T \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/map450_850.R \
        --datafile ${outdir}rds/Noob_${methvalue}_${id}.rds \
        --id after_${methvalue}_${id} \
        --covar F \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

