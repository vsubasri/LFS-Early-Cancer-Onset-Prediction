#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.model.err
#PBS -o log.model.out
#PBS -l vmem=24g
#PBS -l mem=24g

module load R/3.6.1

# Set R library paths - update to your local R installation as needed
# export R_LIBS_USER=/path/to/your/R/library
# export R_LIBS=/path/to/your/R/library

scripts_dir="$(dirname "$0")"

echo "Rscript ${scripts_dir}/sample_refit.R --id $id --aggregate TSS200 --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap --systreat $covars --probes $probes --genes $genes"

Rscript "${scripts_dir}/sample_refit.R" --id $id --aggregate TSS200 --systreat --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit.R" --id $id --aggregate Body --systreat --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit.R" --id $id --aggregate 3UTR --systreat --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit.R" --id $id --aggregate 5UTR --systreat --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit.R" --id $id --aggregate TSS1500 --systreat --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit.R" --id $id --aggregate 1stExon --systreat --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap $covars --probes $probes --genes $genes

