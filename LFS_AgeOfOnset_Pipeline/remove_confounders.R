#####################################################################
############################# Libraries #############################
#####################################################################

require(minfi)
require(ggplot2)
require(dplyr)
require(reshape2)
require(umap)
require(argparse)

source('Scripts/utils.R')

parser <- ArgumentParser()
parser$add_argument("-o","--outdir", action="store")
parser$add_argument("-i","--id", action="store")
parser$add_argument("-d","--data", action="store")

datafile <- args$data
outdir <- args$outdir
id <- args$id 

#####################################################################
############################ Main Script ############################
#####################################################################

## read in data ## 
cat("[ Input data  ]","\n")
data <- read.csv(datafile)

## pca before correction ## 
cat("[ PCA of all data before correction ]","\n")
all_probes <- colnames(data)[grepl('cg',colnames(data))]
pc <- prcomp(as.matrix(data[all_probes]), scale = TRUE)
pc_clin <- cbind(data[,-all_probes],pcx)
write.csv(pc_clin,paste0(outdir,'Output/Noob_',id,'_PCA.csv'),quote=F,row.names=F)
generate_pcsummary(pc_clin,paste0("Noob_",id,"_PCA_summary.csv"),outdir)
generate_pcplots(pc_clin,paste0("Noob_",id),outdir)

## remove outliers ## 
keep <- remove_outliers(pc_clin,3)
data <- data[data$SentrixID %in% keep,]

## array confounder correction ## 
data_450k <- data[data$array == "450",]
if (dim(data_450k) == NULL) {
	cat("No array correction required, skipping to batch correction")
} else {
	## correct 450k data ##
	corrected_450k <- remove_array_confounder(data_450k)
	data_850k <- data[data$array == "850",]
	data <- rbind(corrected_450k,data_850k)
}

## batch confounder correction ## 
data <- remove_batch_confounder(data)
saveRDS(data,paste0(outdir,"rds/NoobCorrected_",id,".rds"))

## pca after correction ##
cat("[ PCA of all data after correction ]","\n")
pc <- prcomp(as.matrix(data[all_probes]), scale = TRUE)
pc_clin <- cbind(data[,-all_probes],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/NoobCorrected_',id,'_PCA.csv'),quote=F,row.names=F)
generate_pcsummary(pc_clin,paste0("NoobCorrected_",id,"_PCA_summary.csv"),outdir)
generate_pcplots(pc_clin,paste0("NoobCorrected_",id),outdir)

