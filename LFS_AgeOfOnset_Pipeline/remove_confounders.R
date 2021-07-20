#####################################################################
############################# Libraries #############################
#####################################################################

require(minfi)
require(ggplot2)
require(dplyr)
require(reshape2)
require(umap)
require(argparse)

source('Scripts/generalUtils.R')

parser <- ArgumentParser()
parser$add_argument("-o","--outdir", action="store")
parser$add_argument("-i","--id", action="store")
parser$add_argument("-d","--data", action="store")

#####################################################################
############################ Main Script ############################
#####################################################################

cat("[ Input data  ]","\n")
data <- read.csv(args$data)

## Predict projection onto PC space in test data

cat("[ Predict PCA in Test Data ]","\n")
data_adj_test <- remove_batch_effect(data,outdir,id)

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


