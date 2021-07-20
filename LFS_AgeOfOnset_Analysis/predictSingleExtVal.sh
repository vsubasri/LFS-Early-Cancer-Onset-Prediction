#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.model.err
#PBS -o log.model.out
#PBS -l vmem=28g
#PBS -l mem=28g
#PBS -l walltime=24:00:00

module load R/3.6.1

export R_LIBS_USER=/hpf/largeprojects/adam/valli/R/3.6.1/
export R_LIBS=/hpf/largeprojects/adam/valli/R/3.6.1/

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/predictSingleExtVal.R --id $id $vars --outdir $outdir

