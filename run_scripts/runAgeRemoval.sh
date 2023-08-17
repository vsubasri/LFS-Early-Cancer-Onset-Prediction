#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time 48:00:00
#SBATCH --job-name ageRemoval

module unload R ; module load R/3.6.1

scripts_dir=/hpf/largeprojects/davidm/vsubasri/methyl_data/Scripts/LFS_ageofonset

Rscript ${scripts_dir}/ageRemoval.R --id NoobCorrected_after_covs_beta_fSVAageofonset_TrainingSet

Rscript ${scripts_dir}/ageRemoval.R --id NoobCorrected_after_covs_beta_fSVAageofonset_TestSet

Rscript ${scripts_dir}/ageRemoval.R --id NoobCorrected_after_covs_beta_fSVAageofonset_ValidationSet

Rscript ${scripts_dir}/ageRemoval.R --id NoobCorrected_after_covs_beta_baPredComBat_TrainingSet

Rscript ${scripts_dir}/ageRemoval.R --id NoobCorrected_after_covs_beta_baPredComBat_TestSet

Rscript ${scripts_dir}/ageRemoval.R --id NoobCorrected_after_covs_beta_baPredComBat_ValidationSet

