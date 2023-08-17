suppressMessages(require('bapred'))
suppressMessages(require('dplyr'))
suppressMessages(require('argparse'))
suppressMessages(require('ggplot2'))
suppressMessages(require('reshape2'))

setwd('/hpf/largeprojects/davidm/vsubasri/methyl_data')
source('Scripts/LFS_ageofonset/util_functions.R')

parser <- ArgumentParser()
parser$add_argument("--id", action="store")
parser$add_argument("--seed", action="store")
parser$add_argument("--outdir", action="store")
parser$add_argument("--nsplit", action="store",type="integer")

args <- parser$parse_args()
id <- args$id
seed <- args$seed
nsplit <- args$nsplit
outdir <- args$outdir

set.seed(seed)


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

cat("[ Run ComBat on training data ]","\n")
Xsubtrain <- as.matrix(data_train[(ind+1):length(data_train)])
batchsubtrain <- factor(as.numeric(factor(data_train$Project,levels=c("SUB14749","SUB14750", "SUB14751","TORONTO","MONTREAL","SUB12648"))),levels=c(1,2,3,4,5,6))
## Save combat parameters from training data
combatparams <- ba(x=Xsubtrain, batch=batchsubtrain, method = "combat")
saveRDS(combatparams,paste0(outdir,'rds/',id,'_baPredComBat_params.rds'))
## Get adjusted data
Xsubtraincombat <- combatparams$xadj
data_Xsubtraincombat <- cbind(data_train[1:ind],Xsubtraincombat)
saveRDS(data_Xsubtraincombat,paste0(outdir,'rds/',id, '_baPredComBat_TrainingSet.rds'))
# Evaluate success of batch effect adjustment on training data
p_train <- prcomp(Xsubtraincombat,scale=T)
p_train <- cbind(data_train[1:ind],p_train$x)
write.csv(p_train,paste0(outdir,'Output/',id, '_baPredComBat_TrainingSet_PCA.csv'))

cat("[ Run ComBat on test data ]","\n")
## Apply correction to test data using training data parameters
Xsubtest <- as.matrix(data_test[(ind+1):length(data_test)])
# Test data before batch effect adjustment 
batchsubtest <- factor(as.numeric(factor(data_test$Project,levels=c("SUB14749","SUB14750", "SUB14751","TORONTO","MONTREAL","SUB12648"))),levels=c(1,2,3,4,5,6))
# Addon batch effect adjustment on test set
Xsubtestcombat <- baaddon(params=combatparams, x=Xsubtest,batch=batchsubtest)
data_Xsubtestcombat <- cbind(data_test[1:ind],Xsubtestcombat)
saveRDS(data_Xsubtestcombat,paste0(outdir,'rds/',id, '_baPredComBat_TestSet.rds'))

cat("[ Run ComBat on validation data ]","\n")
## Apply correction to validation data using training data parameters
Xsubval <- as.matrix(data_val[(ind+1):length(data_val)])
# Validation data before batch effect adjustment
batchsubval <- factor(as.numeric(factor(data_val$Project,levels=c("SUB14749","SUB14750", "SUB14751","TORONTO","MONTREAL","SUB12648"))),levels=c(1,2,3,4,5,6))
# Addon batch effect adjustment on validation set
Xsubvalcombat <- baaddon(params=combatparams, x=Xsubval,batch=batchsubval)
data_Xsubvalcombat <- cbind(data_val[1:ind],Xsubvalcombat)
saveRDS(data_Xsubvalcombat,paste0(outdir, 'rds/',id, '_baPredComBat_ValidationSet.rds'))

# Evaluate success of batch effect adjustment on test data
p_test <- prcomp(Xsubtestcombat,scale=T)
p_test <- cbind(data_test[1:ind],p_test$x)
write.csv(p_test,paste0(outdir,'Output/',id, '_baPredComBat_TestSet_PCA.csv'))

cat("[ Combine and visualize data ]","\n")
# Combine data and plot
allCleaned <- rbind(data_Xsubtraincombat,data_Xsubtestcombat, data_Xsubvalcombat)
saveRDS(allCleaned, paste0(outdir, 'rds/',id, '_baPredComBat.rds'))
pc <- prcomp(allCleaned[(ind+1):length(allCleaned)], scale=T)
pc_clin <- cbind(allCleaned[1:ind],pc$x)
write.csv(p_test,paste0(outdir,'Output/',id, '_baPredComBat_PCA.csv'),quote=F)
generate_pcsummary(pc_clin,paste0(id,'_baPredComBat_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0(id,'_baPredComBat'),outdir)

# Plot technical replicates
duplicated_data <- get_technicalreplicates(allCleaned)
plot_concordance(duplicated_data,paste0(id,"_baPredComBat"),outdir)

pc_beta <- data.frame(duplicated_data[(ind+1):length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir, 'Output/',id,'_baPredComBat_TechnicalReplicates_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0(id,'_baPredComBat_TechnicalReplicates_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0(id,'_baPredComBat_TechnicalReplicates'),outdir)

