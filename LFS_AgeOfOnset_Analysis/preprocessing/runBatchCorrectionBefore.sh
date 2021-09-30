#!/bin/bash
#PBS -S /bin/bash
#PBS -e log/log.BatchCorrection.err
#PBS -o log/log.BatchCorrection.out
#PBS -l vmem=72g
#PBS -l mem=72g
#PBS -l walltime=72:00:00

module load R/3.4.4

export R_LD_LIBRARY_PATH=/hpf/tools/centos6/libgsl/2.1/lib
export R_LIBS_USER=/hpf/largeprojects/adam/valli/R/3.4.4/
export R_LIBS=/hpf/largeprojects/adam/valli/R/3.4.4/

cd /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/map450_850.R \
        --datafile ${outdir}rds/Noob_${methvalue}.rds \
        --id covs_${methvalue} \
        --covar T \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/map450_850.R \
        --datafile ${outdir}rds/Noob_${methvalue}.rds \
        --id ${methvalue} \
        --covar F \
        --seed $seed \
        --outdir ${outdir} \
	--nsplit $nsplit

