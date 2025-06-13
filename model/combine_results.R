library(dplyr)
library(stringr)
library(ggplot2)
library(ROCR)
library(data.table)
library(caret)
library(gridExtra)
library(grid)
library(argparse)

# Use relative paths
script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, 'modelUtils.R'))

parser <- ArgumentParser()
parser$add_argument("--id", action="store")
parser$add_argument("--dir", action="store")
args <- parser$parse_args()
id <- args$id
dir <- args$dir

auc_all_test <- get_all_auc(paste0(dir,"rds/"),"ageofonset_ROCInfoTest.rds")
auc_all_test$dataset <- "test"
auc_all_val <- get_all_auc(paste0(dir,"rds/"),"ageofonset_ROCInfoVal.rds")
auc_all_val$dataset <- "validation"

auc_all <- rbind(auc_all_val, auc_all_test)
write.csv(auc_all, paste0(dir,"Output/",id,'_auc_all.csv'))

