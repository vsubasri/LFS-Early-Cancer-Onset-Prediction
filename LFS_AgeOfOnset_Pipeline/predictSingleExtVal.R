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

outdir <- args$outdir

#####################################################################
############################ Main Script ############################
#####################################################################

## Read in data ## 
data_test <- readRDS(paste0(outdir,'rds/',id,'TestSet.rds'))
## Read in features ## 
features <- read.csv('features.txt',sep='\t')
## Aggregate probes  ## 
data_test <- aggregate_probes(data_test,features)
## Scale data ## 
data_test <- scale_df(data_test,features$gene)

## Predict on new data ## 
xgboost_results <- pred_cancer_xgboost_test(data_test,features$gene)
## Calibrate results ## 
calibrated_results <- platt_scaling(xgboost_results)
## Generate prediction metrics ## 
ROCobj_test <- ROCInfo_atcutoff(calibrated_results,id)
## Save results object ## 
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'_AgeOfOnsetResults.rds'))
