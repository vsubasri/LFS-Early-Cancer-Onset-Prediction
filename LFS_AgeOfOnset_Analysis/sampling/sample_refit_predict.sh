#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.model.err
#PBS -o log.model.out
#PBS -l vmem=8g
#PBS -l mem=8g

module load R/3.6.1

export R_LIBS_USER=/hpf/largeprojects/adam/valli/R/3.6.1/
export R_LIBS=/hpf/largeprojects/adam/valli/R/3.6.1/

echo "Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit_predict.R --id $id --aggregate TSS200 --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap --systreat $covars --probes $probes --genes $genes"

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit_predict.R --id $id --aggregate TSS200 --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit_predict.R --id $id --aggregate Body --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit_predict.R --id $id --aggregate TSS1500 --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit_predict.R --id $id --aggregate 1stExon --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit_predict.R --id $id --aggregate 3UTR --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript /hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/sampling/sample_refit_predict.R --id $id --aggregate 5UTR --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

