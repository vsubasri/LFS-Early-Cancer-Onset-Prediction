suppressMessages(library(minfi))
suppressMessages(library(ggplot2))
suppressMessages(library(gplots))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(umap))
suppressMessages(library(argparse))

# Working directory should be project root when called from shell scripts

parser <- ArgumentParser()
parser$add_argument("--id", action="store")
args <- parser$parse_args()
id <- args$id

set.seed(123)

source('Scripts/util_functions.R')

################################ READ IN INPUT ID FILES #################################

cat("[ Reading in data ]","\n")
data <- readRDS(paste0('rds/',id,'.rds'))
data <- data[!is.na(data$agesamplecollection),]
meth <- data[45:length(data)]

##################################################################################################################
## Regress out cell types ##
###################################################################################################################

if(!file.exists(paste0('rds/',id,'_AgeAdj.rds'))){
	# create a list to store residuals
	resid <- list()
	# loop through all methylation columns and regress each column on the age of sample collection
	for (i in 1:ncol(meth)){

  		temp <- meth[, i]

	# linear model - store residuals
		resid[[i]] <- lm(temp ~ data$agesamplecollection)$residuals

		if (i %% 10000 == 0) {
			print(i)
		}
	}

	adj.meth <- data.frame(do.call('cbind', resid), stringsAsFactors = FALSE)
	colnames(adj.meth) <- colnames(meth)
	dataAgeAdj <- cbind(data[1:44],adj.meth)
	saveRDS(dataAgeAdj,paste0('rds/',id,'_AgeAdj.rds')) 
} else {
	dataAgeAdj <- readRDS(paste0('rds/',id,'_AgeAdj.rds'))
	adj.meth <- dataAgeAdj[45:length(dataAgeAdj)]
}


## generate train/test splits
#data_val <- dataAgeAdj[dataAgeAdj$ids %in% seq(6056,6147,1),]
#saveRDS(data_val,paste0("rds/",id,"_AgeAdj_ValidationSet.rds"))
## create validation set
#dataAgeAdj <- dataAgeAdj[!dataAgeAdj$ids %in% seq(6056,6147,1),]
#testids <- dataAgeAdj[1:44] %>%
#        filter(!ids %in% seq(6056,6147,1)) %>%
#        group_by(ageoutcome = ifelse(dataAgeAdj$ageofonset > 12*6 | is.na(dataAgeAdj$ageofonset),"GT6","LT6")) %>%
#        sample_n(25) %>%
#        select(tm_donor)
#data_test <- dataAgeAdj[dataAgeAdj$tm_donor %in% testids$tm_donor,]
#saveRDS(data_test,paste0("rds/",id,"_AgeAdj_TestSet.rds"))
## create training set
#data_train <- dataAgeAdj[!dataAgeAdj$tm_donor %in% c(data_test$tm_donor,data_val$tm_donor),]

## remove outliers in training data
p_train <- prcomp(data_train[45:length(data_train)],scale=T)
p_train <- cbind(data_train[1:44],p_train$x)
keep_train <- remove_outliers(p_train,3)
data_train <- data_train[data_train$SentrixID %in% keep_train,]
saveRDS(data_train,paste0("rds/",id,"_AgeAdj_TrainingSet.rds"))

## PCA after cell type adjustment ##
cat("[ Plotting corrected data ]","\n")
pc <- prcomp(as.matrix(adj.meth))
pc_clin <- cbind(dataAgeAdj[1:44],pc$x)
write.csv(pc_clin,paste0('Output/',id,'_AgeAdj_PCA.csv'),quote=F)

generate_pcsummary(pc_clin,paste0(id,'_AgeAdj_PCA_summary.csv'))
generate_pcplots(pc_clin,paste0(id,'_AgeAdj'))

cat("[ Plotting technical replicates ]","\n")

duplicated_data <- get_technicalreplicates(dataAgeAdj)
plot_concordance(duplicated_data,paste0(id,'_AgeAdj'))

pc_beta <- data.frame(duplicated_data[45:length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:45],pc$x)
write.csv(pc_clin,paste0('Output/',id,'_AgeAdj_TechnicalReplicates_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0('Output/',id,'_AgeAdj_TechnicalReplicates_PCA_summary.csv'))
generate_pcplots(pc_clin,paste0('Output/',id,'_AgeAdj_TechnicalReplicates'))

u <- umap(adj.meth) ; ud <- cbind(dataAgeAdj[1:44],u$layout)
write.csv(ud,paste0('Output/Umap_',id,'_AgeAdj.csv'))
