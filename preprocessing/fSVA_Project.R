## load up the packages we will need:

suppressMessages(require(argparse))
suppressMessages(require(dplyr))
suppressMessages(require(sva))
suppressMessages(require(ggplot2))
suppressMessages(require(reshape2))

## set up parser

parser <- ArgumentParser()
parser$add_argument("--id", action="store")
parser$add_argument("--outcome", action="store")
parser$add_argument("--seed", action="store")
parser$add_argument("--nsplit", action="store",type="integer")
parser$add_argument("--outdir", action="store")

args <- parser$parse_args()
id <- args$id
outcome <- args$outcome
seed <- args$seed
nsplit <- args$nsplit
outdir <- args$outdir

set.seed(seed)

setwd('/hpf/largeprojects/davidm/vsubasri/methyl_data')
source('Scripts/LFS_ageofonset/util_functions.R')

cat("[ Generate train/test/val split  ]","\n")
data <- readRDS(paste0(outdir,"rds/",id,".rds"))
data <- data[!is.na(data$ids),]
data <- data[!data$Meth %in% c("Problem"),]

datasplits <- getTrainTestVal(data,seed,nsplit)
data_train <- datasplits[[1]]
data_test <- datasplits[[2]]
data_val <- datasplits[[3]]
rm(data)

ind <- grep("cg", colnames(data_train))[1]-1

## Prepare model matrices ##
pheno_train <- data_train[1:ind]
if (outcome == "ageofonset") {
	pheno_train$ageoutcome <- ifelse(pheno_train$ageofonset > 12*6 | is.na(pheno_train$ageofonset),"GT6","LT6")
	trainMod = model.matrix(~ageoutcome,data=pheno_train)
}
if (outcome == "tissue_type") {
	trainMod = model.matrix(~tissue_type,data=pheno_train)
}
if (outcome == "cancer_diagnosis") {
	trainMod = model.matrix(~cancer_diagnosis,data=pheno_train)
}
trainMod0 = model.matrix(~1,data=pheno_train)
trainSv = sva(t(data_train[(ind+1):length(data_train)]),trainMod,trainMod0)

## Run fSVA ##
fsvaobj = fsva(t(data_train[(ind+1):length(data_train)]),trainMod,trainSv,t(data_test[(ind+1):length(data_test)]))
## Format and plot adjusted training data ##
trainCleaned = cbind(data_train[1:ind],t(fsvaobj$db))
saveRDS(trainCleaned,paste0(outdir,'rds/',id,'_fSVA',outcome,'TrainingSet.rds'))
p_train <- prcomp(t(fsvaobj$db))
p_train <- cbind(data_train[1:ind],p_train$x)
write.csv(p_train, paste0(outdir,'Output/',id,'_fSVA',outcome,'TrainingSet_PCA.csv'))

## Format and plot adjusted test data ##
testCleaned = cbind(data_test[1:ind],t(fsvaobj$new))
saveRDS(testCleaned,paste0(outdir,'rds/',id,'_fSVA',outcome,'TestSet.rds'))
p_test <- prcomp(t(fsvaobj$new))
p_test <- cbind(data_test[1:ind],p_test$x)
write.csv(p_test, paste0(outdir,'Output/',id,'_fSVA',outcome,'TestSet_PCA.csv'))

## Save surrogate variables ## 
saveRDS(fsvaobj$newsv,paste0(outdir,'rds/',id,'_fSVA',outcome,'SurrogateVariables.rds'))

## Apply fSVA to external validation set
fsvaobj_val = fsva(t(data_train[(ind+1):length(data_train)]),trainMod,trainSv,t(data_val[(ind+1):length(data_val)]))
valCleaned = cbind(data_val[1:ind],t(fsvaobj_val$new))
saveRDS(valCleaned,paste0(outdir,'rds/',id,'_fSVA',outcome,'ValidationSet.rds'))

## Combine data ##
allCleaned <- rbind(trainCleaned,testCleaned,valCleaned)
saveRDS(allCleaned,paste0(outdir,'rds/',id,'_fSVA',outcome,'.rds'))
pc <- prcomp(allCleaned[(ind+1):length(allCleaned)], scale=T)
pc_clin <- cbind(allCleaned[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/',id,'_fSVA',outcome,'_PCA.csv'),quote=F,row.names=F)
generate_pcsummary(pc_clin,paste0(id,'_fSVA',outcome,'_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0(id,'_fSVA',outcome),outdir)

## Plot concordance between technical replicates ##
duplicated_data <- get_technicalreplicates(allCleaned)
plot_concordance(duplicated_data,paste0(id,"_fSVA",outcome),outdir)

pc_beta <- data.frame(duplicated_data[(ind+1):length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/',id,'_fSVA',outcome,'_TechnicalReplicates_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0(id,'_fSVA',outcome,'_TechnicalReplicates_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0(id,'_fSVA',outcome,'_TechnicalReplicates'),outdir)
