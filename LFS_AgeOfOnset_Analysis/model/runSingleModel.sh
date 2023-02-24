#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.model.err
#PBS -o log.model.out
#PBS -l vmem=48g
#PBS -l mem=48g
#PBS -l walltime=24:00:00

module load R/3.6.1

export R_LIBS_USER=/hpf/largeprojects/adam/valli/R/3.6.1/
export R_LIBS=/hpf/largeprojects/adam/valli/R/3.6.1/


echo "Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.R --id $id $vars --outdir $outdir"

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/runSingleModel.R --id $id $vars --outdir $outdir

