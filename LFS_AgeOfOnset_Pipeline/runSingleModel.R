#####################################################################
############################# Libraries #############################
#####################################################################

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

source('Scripts/modelUtils.R')
set.seed(123)

parser <- ArgumentParser()
parser$add_argument("-o","--outdir", action="store")
parser$add_argument("-i","--id", action="store")

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

data_test <- readRDS(paste0(outdir,'rds/',id,'TestSet.rds'))
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
id <- paste0(id, "scaled_")

id <- paste0(id, "canceratdraw_")
id <- paste0(id, "systreat_")

genes <- colnames(data_train)[45:length(data_train)]

## xgboost
xgboost_results <- pred_cancer_xgboost(data_train,data_test,72,args$sex,args$array,args$canceratdraw,args$systreat,genes)
saveRDS(xgboost_results,paste0(outdir,'rds/',id,'xgboost_ageofonset_results.rds'))
