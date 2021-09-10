#####################################################################
############################# Libraries #############################
#####################################################################

require(stringr)
require(dplyr)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
require(pROC)
require(doParallel)
require(e1071)
require(glmnet)
require(PRROC)
require(ROCR)
require(randomForest)
require(argparse)
require(caret)

source('Scripts/utils.R')
set.seed(123)

## set up parser ##
parser <- ArgumentParser()
parser$add_argument("-o","--outdir", action="store")
parser$add_argument("-i","--id", action="store")

id <- args$id
outdir <- args$outdir

#####################################################################
############################ Main Script ############################
#####################################################################

## Read in data ## 
data_test <- readRDS(paste0(outdir,'rds/',id,'TestSet.rds'))
## Read in model features ## 
probes <- scan('model_probes.txt',what="character",sep='\n')
## Aggregate probes  ## 
data_test <- aggregate_probes(data_test,probes)
## Scale data ## 
data_test <- scale_df(data_test,probes)
id <- paste0(id, "lfs_5UTR_scaled_")

## Read in predictive model ##
xgboost_model <- readRDS('NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_canceratdraw_systreat_xgboost_ageofonset_results.rds')[[1]]
## Predict on new data ## 
xgboost_results <- pred_cancer_xgboost_test(xgboost_model, data_test,72,TRUE,FALSE,TRUE,TRUE,probes)
## Read in platt scaling recalibration model ## 
recal_model <- readRDS('NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_canceratdraw_systreat_xgboost_recalibration_model.rds')
## Calibrate results ## 
calibrated_results <- platt_scaling(recal_model,xgboost_results)
## Generate prediction metrics ## 
ROCobj_test <- ROCInfo_atcutoff(calibrated_results,"test_pred_calibrated","test_label",0.5,paste0(id,"xgboost"),"xgboost")
## Save results object ## 
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'xgboost_ageofonset_ROCInfoTest.rds'))
