setwd('/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data')

#####################################################################
############################# Libraries #############################
#####################################################################

require(stringr)
require(dplyr, warn.conflicts = FALSE)
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
parser$add_argument("-a","--aggregate", action="store")
parser$add_argument("-r", "--array", action="store_true", default=FALSE, help="Include array as a covariate")
parser$add_argument("-x", "--sex", action="store_true", default=FALSE, help="Include sex as a covariate")
parser$add_argument("-c", "--cancer", action="store_true", default=FALSE, help="Remove cancer associated probes")
parser$add_argument("-l", "--lfs", action="store_true", default=FALSE, help="Subset LFS probes")
parser$add_argument("-g", "--gene", action="store_true", default=FALSE, help="Subset probes in gene")
parser$add_argument("-s", "--scale", action="store_true", default=FALSE, help="Scale data")
parser$add_argument("-d", "--canceratdraw", action="store_true", default=FALSE, help="Include cancer at draw variable")
parser$add_argument("-t", "--systreat", action="store_true", default=FALSE, help="Include systemic treatment variable")
args <- parser$parse_args()

id <- args$id
outdir <- args$outdir

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

if (args$lfs) {
  cat(paste0("[ Get LFS probes ]","\n"))
  data_train <- get_lfs_probes(data_train)
  data_test <- get_lfs_probes(data_test)
  id <- paste0(id,"lfs_")
}

if (args$gene) {
  cat(paste0("[ Get probes in gene body ]","\n"))
  data_train <- get_gene_probes(data_train)
  data_test <- get_gene_probes(data_test)
  id <- paste0(id,"gene_")
}

if (args$cancer) {
  cancer_probes <- get_cancer_probes(data_train,id, outdir)
  keep_probes <- colnames(data_train)[!colnames(data_train) %in% cancer_probes]
  if (length(cancer_probes) != 0) {
        cat(paste0("[ Remove cancer signal ]","\n"))
        data_train <- data_train[keep_probes]
        data_test <- data_test[keep_probes]
        id <- paste0(id,"nocancer_")
  } else {
        cat(paste0("[ No significant cancer signal ]","\n"))
        id <- paste0(id,"nocancerfailed_")
  }
}

if (length(args$aggregate)!=0) {
  cat(paste0("[ Aggregate by functional region ]","\n"))
  data_train <- get_func_probes(data_train,args$aggregate)
  data_test <- get_func_probes(data_test,args$aggregate)
  id <- paste0(id,paste0(args$aggregate,"_"))
}

if(args$scale) {
  data_train <- scale_df(data_train)
  data_test <- scale_df(data_test)
  id <- paste0(id, "scaled_")
}

if (args$canceratdraw) {
  id <- paste0(id, "canceratdraw_")
}

if (args$systreat) {
  id <- paste0(id, "systreat_")
}

genes <- colnames(data_train)[45:length(data_train)]

## xgBoost
cat(paste0("[ Predict using xgBoost ]","\n"))
xgboost_model <- readRDS(paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))[[1]]
xgboost_model_val <- readRDS(paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))[[2]]
ROCobj_test <- ROCInfo_atopt(xgboost_model_val,"test_pred","test_label",1, 1,paste0(id,"xgboost"),"xgboost")
cutoff <- ROCobj_test[[5]]
xgboost_results <- pred_cancer_xgboost_test(xgboost_model, data_test,72,TRUE,FALSE,args$canceratdraw,args$systreat,genes)
cat(paste0("[ Get prediction metrics for gbm ]","\n"))
ROCobj_val <- ROCInfo_atcutoff(xgboost_results,"test_pred","test_label",cutoff,paste0(id,"xgboost"),"xgboost")
pr <- platt_scaling(xgboost_model_val,xgboost_results,"xgboost")
test_results <- pr[[1]] ; val_results <- pr[[2]]
ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",0.5,paste0(id,"xgboost"),"xgboost")
ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",0.5,paste0(id,"xgboost"),"xgboost")
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'xgboost_ageofonset_ROCInfoTest.rds'))
saveRDS(ROCobj_val,paste0(outdir,'rds/',id,'xgboost_ageofonset_ROCInfoVal.rds'))

## RANDOM FOREST
cat(paste0("[ Predict using random forest ]","\n"))
rf_model <- readRDS(paste0(outdir,'rds/',id,'rf_ageofonset_results.rds'))[[1]]
rf_model_val <- readRDS(paste0(outdir,'rds/',id,'rf_ageofonset_results.rds'))[[2]]
ROCobj_test <- ROCInfo_atopt(rf_model_val,"test_pred","test_label",1, 1,paste0(id,"rf"),"rf")
cutoff <- ROCobj_test[[5]]
rf_results <- pred_cancer_rf_test(rf_model, data_test,72,TRUE,FALSE,args$canceratdraw,args$systreat,genes)
cat(paste0("[ Get prediction metrics for random forest ]","\n"))
ROCobj_val <- ROCInfo_atcutoff(rf_results,"test_pred","test_label",cutoff,paste0(id,"rf"),"rf")
pr <- platt_scaling(rf_model_val,rf_results,"rf")
test_results <- pr[[1]] ; val_results <- pr[[2]]
ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",0.5,paste0(id,"rf"),"rf")
ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",0.5,paste0(id,"rf"),"rf")
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'rf_ageofonset_ROCInfoTest.rds'))
saveRDS(ROCobj_val,paste0(outdir,'rds/',id,'rf_ageofonset_ROCInfoVal.rds'))

## ENET
cat(paste0("[ Predict using elastic net ]","\n"))
enet_model <- readRDS(paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))[[1]]
enet_model_val <- readRDS(paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))[[2]]
lambda_index <- readRDS(paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))[[4]]
# get prediction metrics for internal test set at optimal cutoff
ROCobj_test <- ROCInfo_atopt(enet_model_val,"test_pred","test_label",1, 1,paste0(id,"enet"),"enet")
cutoff <- ROCobj_test[[5]]
# predict on external validation set
enet_results <- pred_cancer_enet_test(enet_model, lambda_index,data_test,72,TRUE,FALSE,args$canceratdraw,args$systreat,genes)
cat(paste0("[ Get prediction metrics for elastic net ]","\n"))
 get prediction metrics for external validation set aat previously determined cutoff
ROCobj_val <- ROCInfo_atcutoff(enet_results,"test_pred","test_label",cutoff,paste0(id,"enet"),"enet")
pr <- platt_scaling(enet_model_val,enet_results,"enet")
test_results <- pr[[1]] ; val_results <- pr[[2]]
ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",0.5,paste0(id,"enet"),"xgboost")
ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",0.5,paste0(id,"enet"),"xgboost")
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'enet_ageofonset_ROCInfoTest.rds'))
saveRDS(ROCobj_val,paste0(outdir,'rds/',id,'enet_ageofonset_ROCInfoVal.rds'))


## SVM LINEAR
cat(paste0("[ Predict using SVM with linear kernel ]","\n"))
svmlinear_model <- readRDS(paste0(outdir,'rds/',id,'svmLinear_ageofonset_results.rds'))[[1]]
svmlinear_model_val <- readRDS(paste0(outdir,'rds/',id,'svmLinear_ageofonset_results.rds'))[[2]]
ROCobj_test <- ROCInfo_atopt(svmlinear_model_val,"test_pred","test_label",1, 1,paste0(id,"svmLinear"),"svmLinear")
cutoff <- ROCobj_test[[5]]
svmlinear_results <- pred_cancer_svm_test(svmlinear_model, data_test,72,TRUE,FALSE,args$canceratdraw,args$systreat,genes,"svmLinear")
cat(paste0("[ Get prediction metrics for SVM with linear kernel ]","\n"))
ROCobj_val <- ROCInfo_atcutoff(svmlinear_results,"test_pred","test_label",cutoff,paste0(id,"svmLinear"),"svmLinear")
pr <- platt_scaling(svmlinear_model_val,svmlinear_results,"svmLinear")
test_results <- pr[[1]] ; val_results <- pr[[2]]
ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",0.5,paste0(id,"svmLinear"),"svmLinear")
ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",0.5,paste0(id,"svmLinear"),"svmLinear")
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'svmLinear_ageofonset_ROCInfoTest.rds'))
saveRDS(ROCobj_val,paste0(outdir,'rds/',id,'svmLinear_ageofonset_ROCInfoVal.rds'))

## SVM RADIAL
cat(paste0("[ Predict using SVM with radial kernel ]","\n"))
svmradial_model <- readRDS(paste0(outdir,'rds/',id,'svmRadial_ageofonset_results.rds'))[[1]]
svmradial_model_val <- readRDS(paste0(outdir,'rds/',id,'svmRadial_ageofonset_results.rds'))[[2]]
ROCobj_test <- ROCInfo_atopt(svmradial_model_val,"test_pred","test_label",1, 1,paste0(id,"svmRadial"),"svmRadial")
cutoff <- ROCobj_test[[5]]
svmradial_results <- pred_cancer_svm_test(svmradial_model, data_test,72,TRUE,FALSE,args$canceratdraw,args$systreat,genes,"svmRadial")
cat(paste0("[ Get prediction metrics for SVM with radial kernel ]","\n"))
ROCobj_val <- ROCInfo_atcutoff(svmradial_results,"test_pred","test_label",cutoff,paste0(id,"svmRadial"),"svmRadial")
pr <- platt_scaling(svmradial_model_val,svmradial_results,"svmRadial")
test_results <- pr[[1]] ; val_results <- pr[[2]]
ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",0.5,paste0(id,"svmRadial"),"svmRadial")
ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",0.5,paste0(id,"svmRadial"),"svmRadial")
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'svmRadial_ageofonset_ROCInfoTest.rds'))
saveRDS(ROCobj_val,paste0(outdir,'rds/',id,'svmRadial_ageofonset_ROCInfoVal.rds'))

## GBM
cat(paste0("[ Predict using GBM ]","\n"))
gbm_model <- readRDS(paste0(outdir,'rds/',id,'gbm_ageofonset_results.rds'))[[1]]
gbm_model_val <- readRDS(paste0(outdir,'rds/',id,'gbm_ageofonset_results.rds'))[[2]]
ROCobj_test <- ROCInfo_atopt(gbm_model_val,"test_pred","test_label",1, 1,paste0(id,"gbm"),"gbm")
cutoff <- ROCobj_test[[5]]
gbm_results <- pred_cancer_gbm_test(gbm_model, data_test,72,TRUE,FALSE,args$canceratdraw,args$systreat,genes)
cat(paste0("[ Get prediction metrics for gbm ]","\n"))
ROCobj_val <- ROCInfo_atcutoff(gbm_results,"test_pred","test_label",cutoff,paste0(id,"gbm"),"gbm")
pr <- platt_scaling(gbm_model_val,gbm_results,"gbm")
test_results <- pr[[1]] ; val_results <- pr[[2]]
ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",0.5,paste0(id,"gbm"),"gbm")
ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",0.5,paste0(id,"gbm"),"gbm")
saveRDS(ROCobj_test,paste0(outdir,'rds/',id,'gbm_ageofonset_ROCInfoTest.rds'))
saveRDS(ROCobj_val,paste0(outdir,'rds/',id,'gbm_ageofonset_ROCInfoVal.rds'))

