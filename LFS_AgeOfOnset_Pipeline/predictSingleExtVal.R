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

source('Scripts/modelUtils.R')
set.seed(123)

## set up parser ##
parser <- ArgumentParser()
parser$add_argument("-o","--outdir", action="store")
parser$add_argument("-i","--id", action="store")

#####################################################################
############################ Main Script ############################
#####################################################################

cat(paste0("[ Load and clean data ]","\n"))
data_train <- readRDS(paste0(outdir,'rds/',id,'TrainingSet.rds'))
data_train <- remove_age_probes(data_train)
data_train <- remove_xRP(data_train)
data_train <- remove_duplicates(data_train)

data_test <- readRDS(paste0(outdir,'rds/',id,'ValidationSet.rds'))
data_test <- remove_age_probes(data_test)
data_test <- remove_xRP(data_test)
data_test <- remove_duplicates(data_test)

cat(paste0("[ Get LFS probes ]","\n"))
data_train <- get_lfs_probes(data_train)
data_test <- get_lfs_probes(data_test)
id <- paste0(id,"lfs_")

cat(paste0("[ Aggregate by functional region ]","\n"))
data_train <- get_func_probes(data_train,args$aggregate)
data_test <- get_func_probes(data_test,args$aggregate)
id <- paste0(id,paste0(args$aggregate,"_"))


data_train <- scale_df(data_train)
data_test <- scale_df(data_test)
id <- paste0(id, "scaled_canceratdraw_systreat_")

genes <- colnames(data_train)[45:length(data_train)]

## xgBoost
cat(paste0("[ Predict using xgBoost ]","\n"))
xgboost_model <- readRDS(paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))[[1]]
xgboost_model_val <- readRDS(paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))[[2]]
ROCobj_test <- ROCInfo_atopt(xgboost_model_val,"test_pred","test_label",1, 1,paste0(id,"xgboost"),"xgboost")
cutoff <- ROCobj_test[[5]]
xgboost_results <- pred_cancer_xgboost_test(xgboost_model, data_test,72,TRUE,FALSE,TRUE,TRUE,genes)
cat(paste0("[ Get prediction metrics for gbm ]","\n"))
ROCobj_val <- ROCInfo_atcutoff(xgboost_results,"test_pred","test_label",cutoff,paste0(id,"xgboost"),"xgboost")
pr <- platt_scaling(xgboost_model_val,xgboost_results,"xgboost")
test_results <- pr[[1]] ; val_results <- pr[[2]]
ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",0.5,paste0(id,"xgboost"),"xgboost")
ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",0.5,paste0(id,"xgboost"),"xgboost")
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'xgboost_ageofonset_ROCInfoTest.rds'))
saveRDS(ROCobj_val,paste0(outdir,'rds/',id,'xgboost_ageofonset_ROCInfoVal.rds'))
