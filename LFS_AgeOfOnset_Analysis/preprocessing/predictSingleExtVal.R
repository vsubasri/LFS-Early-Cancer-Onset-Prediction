suppressMessages(library(stringr))
suppressMessages(library(dplyr, warn.conflicts = FALSE))
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library(pROC))
suppressMessages(library(doParallel))
suppressMessages(library(e1071))
suppressMessages(library(glmnet))
suppressMessages(library(PRROC))
suppressMessages(library(ROCR))
suppressMessages(library(randomForest))
suppressMessages(library(argparse))
suppressMessages(library(caret))

setwd('/hpf/largeprojects/davidm/vsubasri/methyl_data')
source('Scripts/LFS_ageofonset/modelUtils.R')

parser <- ArgumentParser()
parser$add_argument("-o","--outdir", action="store")
parser$add_argument("-i","--id", action="store")
parser$add_argument("-a","--aggregate", action="store")
parser$add_argument("-e","--age_cutoff", action="store", type="integer", default=72)
parser$add_argument("-m", "--family", action="store_true", default=FALSE, help="Include family variable")
parser$add_argument("-c", "--cancer", action="store_true", default=FALSE, help="Remove cancer associated probes")
parser$add_argument("-l", "--lfs", action="store_true", default=FALSE, help="Subset LFS probes")
parser$add_argument("-g", "--gene", action="store_true", default=FALSE, help="Subset probes in gene")
parser$add_argument("-s", "--scale", action="store_true", default=FALSE, help="Scale data")
parser$add_argument("-x", "--sex",  action="store_true", default=FALSE, help="Include sex variable")
parser$add_argument("-d", "--canceratdraw", action="store_true", default=FALSE, help="Include cancer at draw variable")
parser$add_argument("-t", "--systreat", action="store_true", default=FALSE, help="Include systemic treatment variable")
parser$add_argument("-p", "--cellprops", action="store_true", default=FALSE, help="Include cell proportion variables")
parser$add_argument("-r", "--array", action="store_true", default=FALSE, help="Include array variable")
parser$add_argument("-w", "--fnw", action="store", help="False negative weight", type="integer")
parser$add_argument("-b", "--threshold", action="store", help="Prediction probability threshold")
parser$add_argument("-f", "--tp53", action="store_true", default=FALSE, help="Include TP53 exon")
args <- parser$parse_args()

id <- args$id
outdir <- args$outdir
FN_weight <- args$fnw

if (!is.null(args$threshold)) {
  cutoff <- args$threshold
}

cat(paste0("[ Load and clean data ]","\n"))
data_train <- readRDS(paste0(outdir,'rds/',id,'TrainingSet.rds'))
data_train <- remove_age_probes(data_train)
data_train <- remove_xRP(data_train)
data_train <- remove_duplicates(data_train)

data_test <- readRDS(paste0(outdir,'rds/',id,'TestSet.rds'))
data_test <- remove_age_probes(data_test)
data_test <- remove_xRP(data_test)
data_test <- remove_duplicates(data_test)

ind <- grep("cg", colnames(data_train))[1]

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

if (args$cellprops) {
  data_train <- get_cell_props(data_train)
  data_test <- get_cell_props(data_test)
  id <- paste0(id, "cellprops_")
}

if (args$family) {
  data_train <- get_family_vars(data_train)
  data_test <- get_family_vars(data_test)
  id <- paste0(id, "family_")
  ind <- ind + 2
}

if (args$scale) {
  data_train <- scale_df(data_train)
  data_test <- scale_df(data_test)
  id <- paste0(id, "scaled_")
}

if (args$sex) {
  id <- paste0(id, "sex_")
}

if (args$canceratdraw) {
  id <- paste0(id, "canceratdraw_")
}

if (args$systreat) {
  id <- paste0(id, "systreat_")
}


if (args$array) {
  id <- paste0(id, "array_")
}

if (args$tp53) {
  id <- paste0(id, "tp53_")
  ind <- ind + 9
}

genes <- colnames(data_train)[ind:length(data_train)]
ageofonset <- args$age_cutoff

id <- paste0(id, ageofonset,"_")

## xgBoost
if (file.exists(paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))) {
  cat(paste0("[ Predict using xgBoost ]","\n"))
  xgboost_model <- readRDS(paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))[[1]]
  xgboost_model_val <- readRDS(paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))[[2]]
  xgboost_results_test <- pred_cancer_xgboost_test(xgboost_model, data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes)
  cat(paste0("[ Get calibrated predictions ]","\n"))
  pr <- platt_scaling(xgboost_model_val,xgboost_results_test,"xgboost")
  val_results <- pr[[1]] ; test_results <- pr[[2]]
  cat(paste0("[ Get prediction metrics for xgboost ]","\n"))
  if (!is.null(args$fnw)) {
    fid <- paste0(id, "FNW",FN_weight,"_")
    ROCobj_val <- ROCInfo_atopt(val_results,"test_pred_calibrated","test_label",1,FN_weight,paste0(fid,"xgboost"),"xgboost")
    cutoff <- ROCobj_val[[5]]
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"xgboost"),"xgboost")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'xgboost_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'xgboost_ageofonset_ROCInfoVal.rds'))
  } else if (!is.null(args$threshold)) {
    fid <- paste0(id, "threshold",cutoff,"_")
    ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"xgboost"),"xgboost")
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"xgboost"),"xgboost")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'xgboost_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'xgboost_ageofonset_ROCInfoVal.rds'))
  }
}
## GBM
if (file.exists(paste0(outdir,'rds/',id,'gbm_ageofonset_results.rds'))) {
  cat(paste0("[ Predict using GBM ]","\n"))
  gbm_model <- readRDS(paste0(outdir,'rds/',id,'gbm_ageofonset_results.rds'))[[1]]
  gbm_model_val <- readRDS(paste0(outdir,'rds/',id,'gbm_ageofonset_results.rds'))[[2]]
  gbm_results_test <- pred_cancer_gbm_test(gbm_model, data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes)
  cat(paste0("[ Get prediction metrics for gbm ]","\n"))
  pr <- platt_scaling(gbm_model_val,gbm_results_test,"gbm")
  val_results <- pr[[1]] ; test_results <- pr[[2]]
  if (!is.null(args$fnw)) {
    fid <- paste0(id, "FNW",FN_weight,"_")
    ROCobj_val <- ROCInfo_atopt(val_results,"test_pred_calibrated","test_label",1,FN_weight,paste0(fid,"gbm"),"gbm")
    cutoff <- ROCobj_val[[5]]
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"gbm"),"gbm")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'gbm_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'gbm_ageofonset_ROCInfoVal.rds'))
  } else if (!is.null(args$threshold)) {
    fid <- paste0(id, "threshold",cutoff,"_")
    ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"gbm"),"gbm")
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"gbm"),"gbm")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'gbm_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'gbm_ageofonset_ROCInfoVal.rds'))
  }
}

## RANDOM FOREST
if (file.exists(paste0(outdir,'rds/',id,'rf_ageofonset_results.rds'))) {
  cat(paste0("[ Predict using random forest ]","\n"))
  rf_model <- readRDS(paste0(outdir,'rds/',id,'rf_ageofonset_results.rds'))[[1]]
  rf_model_val <- readRDS(paste0(outdir,'rds/',id,'rf_ageofonset_results.rds'))[[2]]
  rf_results_test <- pred_cancer_rf_test(rf_model, data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes)
  cat(paste0("[ Get prediction metrics for random forest ]","\n"))
  pr <- platt_scaling(rf_model_val,rf_results_test,"rf")
  val_results <- pr[[1]] ; test_results <- pr[[2]]
  if (!is.null(args$fnw)) {
    fid <- paste0(id, "FNW",FN_weight,"_")
    ROCobj_val <- ROCInfo_atopt(val_results,"test_pred_calibrated","test_label",1,FN_weight,paste0(fid,"rf"),"rf")
    cutoff <- ROCobj_val[[5]]
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"rf"),"rf")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'rf_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'rf_ageofonset_ROCInfoVal.rds'))
  } else if (!is.null(args$threshold)) {
    fid <- paste0(id, "threshold",cutoff,"_")
    ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"rf"),"rf")
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"rf"),"rf")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'rf_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'rf_ageofonset_ROCInfoVal.rds'))
  }
}

## ENET
if (file.exists(paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))) {
  cat(paste0("[ Predict using elastic net ]","\n"))
  enet_model <- readRDS(paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))[[1]]
  enet_model_val <- readRDS(paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))[[2]]
  lambda_index <- readRDS(paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))[[4]]
  # get prediction metrics for internal test set at optimal cutoff
  # predict on external validation set
  enet_results_test <- pred_cancer_enet_test(enet_model, lambda_index,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes)
  cat(paste0("[ Get prediction metrics for elastic net ]","\n"))
  # get prediction metrics for external validation set aat previously determined cutoff
  pr <- platt_scaling(enet_model_val,enet_results_test,"enet")
  val_results <- pr[[1]] ; test_results <- pr[[2]]
  if (!is.null(args$fnw)) {
    fid <- paste0(id, "FNW",FN_weight,"_")
    ROCobj_val <- ROCInfo_atopt(val_results,"test_pred_calibrated","test_label",1,FN_weight,paste0(fid,"enet"),"enet")
    cutoff <- ROCobj_val[[5]]
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"enet"),"enet")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'enet_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'enet_ageofonset_ROCInfoVal.rds'))
  } else if (!is.null(args$threshold)) {
    fid <- paste0(id, "threshold",cutoff,"_")
    ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"enet"),"enet")
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"enet"),"enet")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'enet_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'enet_ageofonset_ROCInfoVal.rds'))
  }
}

## SVM LINEAR
if (file.exists(paste0(outdir,'rds/',id,'svmLinear_ageofonset_results.rds'))) {
  cat(paste0("[ Predict using SVM with linear kernel ]","\n"))
  svmlinear_model <- readRDS(paste0(outdir,'rds/',id,'svmLinear_ageofonset_results.rds'))[[1]]
  svmlinear_model_val <- readRDS(paste0(outdir,'rds/',id,'svmLinear_ageofonset_results.rds'))[[2]]
  svmlinear_results_test <- pred_cancer_svm_test(svmlinear_model, data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes,"svmLinear")
  cat(paste0("[ Get prediction metrics for SVM with linear kernel ]","\n"))
  pr <- platt_scaling(svmlinear_model_val,svmlinear_results_test,"svmLinear")
  val_results <- pr[[1]] ; test_results <- pr[[2]]
  if (!is.null(args$fnw)) {
    fid <- paste0(id, "FNW",FN_weight,"_")
    ROCobj_val <- ROCInfo_atopt(val_results,"test_pred_calibrated","test_label",1,FN_weight,paste0(fid,"svmLinear"),"svmLinear")
    cutoff <- ROCobj_val[[5]]
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmLinear"),"svmLinear")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'svmLinear_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'svmLinear_ageofonset_ROCInfoVal.rds'))
  } else if (!is.null(args$threshold)) {
    fid <- paste0(id, "threshold",cutoff,"_")
    ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmLinear"),"svmLinear")
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmLinear"),"svmLinear")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'svmLinear_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'svmLinear_ageofonset_ROCInfoVal.rds'))
  }
}

## SVM RADIAL
if (file.exists(paste0(outdir,'rds/',id,'svmRadial_ageofonset_results.rds'))) {
  cat(paste0("[ Predict using SVM with radial kernel ]","\n"))
  svmradial_model <- readRDS(paste0(outdir,'rds/',id,'svmRadial_ageofonset_results.rds'))[[1]]
  svmradial_model_val <- readRDS(paste0(outdir,'rds/',id,'svmRadial_ageofonset_results.rds'))[[2]]
  svmradial_results_test <- pred_cancer_svm_test(svmradial_model, data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53, genes,"svmRadial")
  cat(paste0("[ Get prediction metrics for SVM with radial kernel ]","\n"))
  pr <- platt_scaling(svmradial_model_val,svmradial_results_test,"svmRadial")
  val_results <- pr[[1]] ; test_results <- pr[[2]]
  if (!is.null(args$fnw)) {
    fid <- paste0(id, "FNW",FN_weight,"_")
    ROCobj_val <- ROCInfo_atopt(val_results,"test_pred_calibrated","test_label",1,FN_weight,paste0(fid,"svmRadial"),"svmRadial")
    cutoff <- ROCobj_val[[5]]
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmRadial"),"svmRadial")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'svmRadial_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'svmRadial_ageofonset_ROCInfoVal.rds'))
  } else if (!is.null(args$threshold)) {
    fid <- paste0(id, "threshold",cutoff,"_")
    ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmRadial"),"svmRadial")
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmRadial"),"svmRadial")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'svmRadial_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'svmRadial_ageofonset_ROCInfoVal.rds'))
  }
}

## NNET
if (file.exists(paste0(outdir,'rds/',id,'nnet_ageofonset_results.rds'))) {
    cat(paste0("[ Predict using neural network ]","\n"))
    nnet_model <- readRDS(paste0(outdir,'rds/',id,'nnet_ageofonset_results.rds'))[[1]]
    nnet_model_val <- readRDS(paste0(outdir,'rds/',id,'nnet_ageofonset_results.rds'))[[2]]
    nnet_results_test <- pred_cancer_nnet_test(nnet_model, data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes)
    cat(paste0("[ Get prediction metrics for neural network ]","\n"))
    pr <- platt_scaling(nnet_model_val,nnet_results_test,"nnet")
    val_results <- pr[[1]] ; test_results <- pr[[2]]
    if (!is.null(args$fnw)) {
      fid <- paste0(id, "FNW",FN_weight,"_")
      ROCobj_val <- ROCInfo_atopt(val_results,"test_pred_calibrated","test_label",1,FN_weight,paste0(fid,"nnet"),"nnet")
      cutoff <- ROCobj_val[[5]]
      ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"nnet"),"nnet")
      saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'nnet_ageofonset_ROCInfoTest.rds'))
      saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'nnet_ageofonset_ROCInfoVal.rds'))
    } else if (!is.null(args$threshold)) {
      fid <- paste0(id, "threshold",cutoff,"_")
      ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"nnet"),"nnet")
      ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"nnet"),"nnet")
      saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'nnet_ageofonset_ROCInfoTest.rds'))
      saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'nnet_ageofonset_ROCInfoVal.rds'))
    }
}

## SVM POLYNOMIAL
if (file.exists(paste0(outdir,'rds/',id,'svmPoly_ageofonset_results.rds'))) {
  cat(paste0("[ Predict using SVM with polynomial kernel ]","\n"))
  svmpoly_model <- readRDS(paste0(outdir,'rds/',id,'svmPoly_ageofonset_results.rds'))[[1]]
  svmpoly_model_val <- readRDS(paste0(outdir,'rds/',id,'svmPoly_ageofonset_results.rds'))[[2]]
  svmpoly_results_test <- pred_cancer_svm_test(svmpoly_model, data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes,"svmPoly")
  cat(paste0("[ Get prediction metrics for SVM with polynomial kernel ]","\n"))
  pr <- platt_scaling(svmpoly_model_val,svmpoly_results_test,"svmPoly")
  val_results <- pr[[1]] ; test_results <- pr[[2]]
  if (!is.null(args$fnw)) {
    fid <- paste0(id, "FNW",FN_weight,"_")
    ROCobj_val <- ROCInfo_atopt(val_results,"test_pred_calibrated","test_label",1,FN_weight,paste0(fid,"svmPoly"),"svmPoly")
    cutoff <- ROCobj_val[[5]]
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmPoly"),"svmPoly")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'svmPoly_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'svmPoly_ageofonset_ROCInfoVal.rds'))
  } else if (!is.null(args$threshold)) {
    fid <- paste0(id, "threshold",cutoff,"_")
    ROCobj_val <- ROCInfo_atcutoff(val_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmPoly"),"svmPoly")
    ROCobj_test <- ROCInfo_atcutoff(test_results,"test_pred_calibrated","test_label",cutoff,paste0(fid,"svmPoly"),"svmPoly")
    saveRDS(ROCobj_test,paste0(outdir,'rds/',fid,'svmPoly_ageofonset_ROCInfoTest.rds'))
    saveRDS(ROCobj_val,paste0(outdir,'rds/',fid,'svmPoly_ageofonset_ROCInfoVal.rds'))
  }
}
