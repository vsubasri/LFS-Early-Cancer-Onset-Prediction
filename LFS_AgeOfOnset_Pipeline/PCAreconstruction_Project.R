setwd('/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/')

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

#####################################################################
############################ Main Script ############################
#####################################################################

ind <- 44

cat("[ Generate train/test/val split  ]","\n")
data <- readRDS(paste0(outdir,"rds/",id,".rds"))
data <- data[!is.na(data$ids),]

p_train <- readRDS('rds/batchEffectRemoval.rds')

## Predict projection onto PC space in test data
cat("[ Predict PCA in Test Data ]","\n")
test_pred <-predict(p_train,data_test[(ind+1):length(data_test)])
Xhat_pred <- test_pred[, !(colnames(test_pred) %in% c(maxpc))] %*% t(p_train$rotation[, !(colnames(p_train$rotation) %in% c(maxpc))])
beta_adj_pred <- scale(Xhat_pred, center = -(colMeans(data_test[(ind+1):length(data_test)])), scale = T)
data_adj_test <- cbind(data_test[1:ind],beta_adj_pred)
saveRDS(data_adj_test,paste0(outdir, 'rds/',id,'_ProjPC2Adj_TestSet.rds'))

## Predict projection onto PC space in validation data
cat("[ Predict PCA in Validation Data ]","\n")
val_pred <-predict(p_train,data_val[(ind+1):length(data_val)])
Xhat_pred_val <- val_pred[, !(colnames(val_pred) %in% c(maxpc))] %*% t(p_train$rotation[, !(colnames(p_train$rotation) %in% c(maxpc))])
beta_adj_pred_val <- scale(Xhat_pred_val, center = -(colMeans(data_val[(ind+1):length(data_val)])), scale = T)
data_adj_val <- cbind(data_val[1:ind],beta_adj_pred_val)
saveRDS(data_adj_val,paste0(outdir,'rds/',id,'_ProjPC2Adj_ValidationSet.rds'))

## Plot test set correction ##
cat("[ Plotting PCA of corrected test data ]","\n")
pc <- prcomp(beta_adj_pred, scale=T)
pc_clin <- cbind(data_adj_test[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/',id,'_ProjPC2Adj_TestSet_PCA.csv'),quote=F)

## Combined Data ##
allCleaned <- rbind(data_adj, data_adj_test,data_adj_val)
saveRDS(allCleaned,paste0(outdir,'rds/',id,'_ProjPC2Adj.rds'))
pc <- prcomp(allCleaned[(ind+1):length(allCleaned)], scale=T)
pc_clin <- cbind(allCleaned[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/',id,'_ProjPC2Adj_PCA.csv'),quote=F,row.names=F)
generate_pcsummary(pc_clin,paste0(id,'_ProjPC2Adj_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0(id,'_ProjPC2Adj'),outdir)

## Plot concordance between technical replicates ##
duplicated_data <- get_technicalreplicates(allCleaned)
plot_concordance(duplicated_data,paste0(id,"_ProjPC2Adj"),outdir)

pc_beta <- data.frame(duplicated_data[(ind+1):length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/',id,'_ProjPC2Adj_TechnicalReplicates_PCA.csv'),quote=F,row.names=F)
generate_pcsummary(pc_clin,paste0(id,'_ProjPC2Adj_TechnicalReplicates_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0(id,'_ProjPC2Adj_TechnicalReplicates'),outdir)


