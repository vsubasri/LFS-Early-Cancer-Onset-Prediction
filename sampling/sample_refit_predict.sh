#!/bin/bash
#PBS -S /bin/bash
#PBS -e log.model.err
#PBS -o log.model.out
#PBS -l vmem=8g
#PBS -l mem=8g

module load R/3.6.1

# Set R library paths - update to your local R installation as needed
# export R_LIBS_USER=/path/to/your/R/library
# export R_LIBS=/path/to/your/R/library

scripts_dir="$(dirname "$0")"

echo "Rscript ${scripts_dir}/sample_refit_predict.R --id $id --aggregate TSS200 --canceratdraw --model $model --outdir $outdir --bootstrap $bootstrap --systreat $covars --probes $probes --genes $genes"

Rscript "${scripts_dir}/sample_refit_predict.R" --id $id --aggregate TSS200 --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit_predict.R" --id $id --aggregate Body --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit_predict.R" --id $id --aggregate TSS1500 --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit_predict.R" --id $id --aggregate 1stExon --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit_predict.R" --id $id --aggregate 3UTR --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

Rscript "${scripts_dir}/sample_refit_predict.R" --id $id --aggregate 5UTR --outdir $outdir --bootstrap $bootstrap --canceratdraw --model $model --systreat $covars --probes $probes --genes $genes

