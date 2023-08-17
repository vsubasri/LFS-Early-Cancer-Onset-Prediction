suppressMessages(library(stringr))
suppressMessages(library(dplyr))
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
parser$add_argument("-f", "--tp53", action="store_true", default=FALSE, help="Include TP53 exon")
args <- parser$parse_args()

id <- args$id
outdir <- args$outdir

cat(paste0("[ Load and clean data ]","\n"))

data_train <- readRDS(paste0(outdir,'rds/',id,'TrainingSet.rds'))
data_train <- remove_age_probes(data_train)
data_train <- remove_xRP(data_train)
data_train <- remove_duplicates(data_train)

data_test <- readRDS(paste0(outdir,'rds/',id,'ValidationSet.rds'))
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
  cancer_probes <- get_cancer_probes(data_train,id,outdir)
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

cat("[ Train model ]\n")
## xgboost
if (!file.exists(paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))) {
  xgboost_results <- pred_cancer_xgboost(data_train,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops, args$family, args$tp53, genes)
  saveRDS(xgboost_results,paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))
}

## gbm 
if (!file.exists(paste0(outdir,'rds/',id,'gbm_ageofonset_results.rds'))) {
  gbm_results <- pred_cancer_gbm(data_train,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family, args$tp53, genes)
  saveRDS(gbm_results,paste0(outdir,'rds/',id,'gbm_ageofonset_results.rds'))
}

if (args$canceratdraw) {
  data_train <- data_train[!is.na(data_train[c("cancer_atdraw")]),]
  data_test <- data_test[!is.na(data_test[c("cancer_atdraw")]),]
}

if (args$systreat) {
  data_train <- data_train[!is.na(data_train[c("systemic_treatment_atdraw")]),]
  data_test <- data_test[!is.na(data_test[c("systemic_treatment_atdraw")]),]
}

## random forest
if (!file.exists(paste0(outdir,'rds/',id,'rf_ageofonset_results.rds'))) {
  rf_results <- pred_cancer_rf(data_train,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53, genes)
  saveRDS(rf_results,paste0(outdir,'rds/',id,'rf_ageofonset_results.rds'))
}
## enet
if (!file.exists(paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))) {
  enet_results <- pred_cancer_enet(data_train,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53, genes)
  saveRDS(enet_results,paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))
}

## svm linear
if (!file.exists(paste0(outdir,'rds/',id,'svmLinear_ageofonset_results.rds'))) {
  svmLinear_results <- pred_cancer_svm(data_train,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family, args$tp53,genes,"svmLinear")
  saveRDS(svmLinear_results,paste0(outdir,'rds/',id,'svmLinear_ageofonset_results.rds'))
}
## svm radial
if (!file.exists(paste0(outdir,'rds/',id,'svmRadial_ageofonset_results.rds'))) {
  svmRadial_results <- pred_cancer_svm(data_train,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes,"svmRadial")
  saveRDS(svmRadial_results,paste0(outdir,'rds/',id,'svmRadial_ageofonset_results.rds'))
}

## svm polynomial
if (!file.exists(paste0(outdir,'rds/',id,'svmPoly_ageofonset_results.rds'))) {
  svmPoly_results <- pred_cancer_svm(data_train,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops,args$family,args$tp53,genes,"svmPoly")
  saveRDS(svmPoly_results,paste0(outdir,'rds/',id,'svmPoly_ageofonset_results.rds'))
}

## nnet
if (!file.exists(paste0(outdir,'rds/',id,'nnet_ageofonset_results.rds'))) {
  nnet_results <- pred_cancer_nnet(data_train,data_test,ageofonset,args$sex,args$array,args$canceratdraw,args$systreat,args$cellprops, args$family,args$tp53, genes)
  saveRDS(nnet_results,paste0(outdir,'rds/',id,'nnet_ageofonset_results.rds'))
}

