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

# Working directory should be project root when called from shell scripts
source('model/modelUtils.R')

parser <- ArgumentParser()
parser$add_argument("-o","--outdir", action="store")
parser$add_argument("-i","--id", action="store")
parser$add_argument("-a","--aggregate", action="store")
parser$add_argument("-r", "--array", action="store_true", default=FALSE, help="Include array as a covariate")
parser$add_argument("-x", "--sex", action="store_true", default=FALSE, help="Include sex as a covariate")
parser$add_argument("-c", "--cancer", action="store_true", default=FALSE, help="Remove cancer associated probes")
parser$add_argument("-g", "--gene", action="store_true", default=FALSE, help="Subset probes in gene")
parser$add_argument("-s", "--scale", action="store_true", default=FALSE, help="Scale data")
parser$add_argument("-d", "--canceratdraw", action="store_true", default=FALSE, help="Include cancer at draw variable")
parser$add_argument("-t", "--systreat", action="store_true", default=FALSE, help="Include systemic treatment at draw variable")
parser$add_argument("-m","--model", action="store")
parser$add_argument("-b","--bootstrap", action="store")
parser$add_argument("-p","--probes", action="store",type="integer")
parser$add_argument("-f","--genes", action="store",type="integer")
args <- parser$parse_args()

id <- args$id
outdir <- args$outdir

set.seed(123)

cat(paste0("[ Load and clean data ]","\n"))
data_train <- readRDS(paste0(outdir,'rds/',id,'TrainingSet.rds'))
data_train <- remove_age_probes(data_train)
data_train <- remove_xRP(data_train)
data_train <- remove_duplicates(data_train)

data_test <- readRDS(paste0(outdir,'rds/',id,'TestSet.rds'))
data_test <- remove_age_probes(data_test)
data_test <- remove_xRP(data_test)
data_test <- remove_duplicates(data_test)

id <- paste0(id,"bs",args$bootstrap,"_")

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
  rm(keep_probes); rm(cancer_probes)
}

if (args$gene) {
  cat(paste0("[ Get probes in gene body ]","\n"))
  data_train <- get_gene_probes(data_train)
  data_test <- get_gene_probes(data_test)
  id <- paste0(id,"gene_")
}

if (length(args$aggregate)!=0) {
  cat(paste0("[ Sampling probes ]","\n"))
  ## sample probes
  data_train <- sample_func_probes(data_train,args$aggregate,args$genes,args$probes,args$bootstrap)
  data_test <- data_test[colnames(data_train)]
  print(dim(data_train))
  ## aggregate probes
  cat(paste0("[ Aggregate by functional region ]","\n"))
  data_train <- aggregate_sampled_probes(data_train)
  data_test <- aggregate_sampled_probes(data_test)
  id <- paste0(id,paste0(args$aggregate,"_",args$probes,"_",args$genes,"_"))
}

if(args$scale) {
  data_train <- scale_df(data_train)
  data_test <- scale_df(data_test)
  id <- paste0(id, "scaled_")
}

if (args$canceratdraw ) {
 id <- paste0(id, "canceratdraw_") 
}

if (args$systreat) {
  id <- paste0(id, "systreat_")
}

genes <- colnames(data_train)[45:length(data_train)]
## enet
if (args$model == "enet") {
  enet_results <- pred_cancer_enet(data_train,data_test,72,args$sex,args$array,args$canceratdraw,args$systreat,genes)
  saveRDS(enet_results,paste0(outdir,'rds/',id,'enet_ageofonset_results.rds'))
}
## random forest
if (args$model == "randomforest") {
  rf_results <- pred_cancer_rf(data_train,data_test,72,args$sex,args$array,args$canceratdraw,args$systreat,genes)
  saveRDS(rf_results,paste0(outdir,'rds/',id,'rf_ageofonset_results.rds'))
}

## svm linear
if (args$model == "svmLinear") {
  svmLinear_results <- pred_cancer_svm(data_train,data_test,72,args$sex,args$array,args$canceratdraw,args$systreat,genes,"svmLinear")
  saveRDS(svmLinear_results,paste0(outdir,'rds/',id,'svmLinear_ageofonset_results.rds'))
}
## svm radial
if (args$model == "svmRadial") {
  svmRadial_results <- pred_cancer_svm(data_train,data_test,72,args$sex,args$array,args$canceratdraw,args$systreat,genes,"svmRadial")
  saveRDS(svmRadial_results,paste0(outdir,'rds/',id,'svmRadial_ageofonset_results.rds'))
}
## gbm 
if (args$model == "gbm") {
  gbm_results <- pred_cancer_gbm(data_train,data_test,72,args$sex,args$array,args$canceratdraw,args$systreat,genes)
  saveRDS(gbm_results,paste0(outdir,'rds/',id,'gbm_ageofonset_results.rds'))
}
## xgboost
if (args$model == "xgboost") {
  xgboost_results <- pred_cancer_xgboost(data_train,data_test,72,args$sex,args$array,args$canceratdraw,args$systreat,genes)
  saveRDS(xgboost_results,paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))
}

