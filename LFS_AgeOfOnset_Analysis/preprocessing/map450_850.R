setwd('/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/')

#####################################################################
############################# Libraries #############################
#####################################################################

require(ggplot2)
require(argparse)
require(dplyr)
require(limma)
require(data.table)
require(umap)
require(tidyr)

## set up parser ##
parser <- ArgumentParser()
parser$add_argument("--datafile", action="store")
parser$add_argument("--id", action="store")
parser$add_argument("--covar", action="store")
parser$add_argument("--seed", action="store")
parser$add_argument("--nsplit", action="store",type="integer")
parser$add_argument("--outdir", action="store")

## set variables from parser ##
args <- parser$parse_args()
datafile <- args$datafile
covar <- args$covar
seed <- args$seed
nsplit <- args$nsplit
outdir <- args$outdir
id <- args$id

source('Scripts/generalUtils.R')
set.seed(seed)

#####################################################################
############################ Main Script ############################
#####################################################################

## read in input data ##
cat("[ Reading in input data ]","\n")
data <- readRDS(datafile)
data <- data[!is.na(data$cancer_diagnosis),]
data$array <- factor(data$array)

## remove outliers ## 
cat("[ PCA of all data before correction ]","\n")
pc <- prcomp(as.matrix(data[45:length(data)]), scale = TRUE)
pc_clin <- cbind(data[1:44],pc$x)
keep <- remove_outliers(pc_clin,3)
data <- data[data$SentrixID %in% keep,]
data_450k <- data[data$array == "450",]

## correct 450k data ##
corrected_450k <- run_correction(data,id,covar,outdir)
data_850k <- data[data$array == "850",]
all_corrected <- rbind(corrected_450k,data_850k)
all_corrected <- all_corrected[!is.na(all_corrected[length(all_corrected)]),]
saveRDS(all_corrected,paste0(outdir,"rds/NoobCorrected_",id,".rds"))

## plot correction ##
cat("[ PCA of all data after correction ]","\n")
pc <- prcomp(as.matrix(all_corrected[45:length(all_corrected)]), scale = TRUE)
pc_clin <- cbind(all_corrected[1:44],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/NoobCorrected_',id,'_PCA.csv'),quote=F,row.names=F)
generate_pcsummary(pc_clin,paste0("NoobCorrected_",id,"_PCA_summary.csv"),outdir)
generate_pcplots(pc_clin,paste0("NoobCorrected_",id),outdir)

u <- umap(all_corrected[45:length(all_corrected)]) ; ud <- cbind(all_corrected[1:44],u$layout)
write.csv(ud,paste0(outdir,'Output/Umap_NoobCorrected_',id,'.csv'))

## plot correction for technical replicates only ##
cat("[ PCA of all technical replicates after correction ]","\n")
duplicated_data <- get_technicalreplicates(all_corrected)
plot_concordance(duplicated_data,paste0("NoobCorrected_",id),outdir)

pc_beta <- data.frame(duplicated_data[45:length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:44],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/NoobCorrected_',id,'_TechnicalReplicates_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0('NoobCorrected_',id,'_TechnicalReplicates_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0('NoobCorrected_',id,'_TechnicalReplicates'),outdir)

## generate train/test/val splits
datasplits <- getTrainTestVal(all_corrected,seed,nsplit)

data_train <- datasplits[[1]]
saveRDS(data_train,paste0(outdir,"rds/NoobCorrected_",id,"_TrainingSet.rds"))
data_test <- datasplits[[2]]
saveRDS(data_test,paste0(outdir,"rds/NoobCorrected_",id,"_TestSet.rds"))
data_val <- datasplits[[3]]
saveRDS(data_val,paste0(outdir,"rds/NoobCorrected_",id,"_ValidationSet.rds"))


