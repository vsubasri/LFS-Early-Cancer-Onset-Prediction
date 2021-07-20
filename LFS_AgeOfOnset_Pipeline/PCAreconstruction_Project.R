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
beta <- data[ , grepl( "cg" , names( data ) ) ]
clin <- data[ , -grepl( "cg" , names( data ) ) ]

p_train <- readRDS('rds/batchEffectRemoval.rds')

## Predict projection onto PC space in test data
cat("[ Predict PCA in Test Data ]","\n")
test_pred <-predict(p_train,beta)
Xhat_pred <- test_pred[, !(colnames(test_pred) %in% c(maxpc))] %*% t(p_train$rotation[, !(colnames(p_train$rotation) %in% c(maxpc))])
beta_adj_pred <- scale(Xhat_pred, center = -(colMeans(beta)), scale = T)
data_adj_test <- cbind(clin,beta_adj_pred)
saveRDS(data_adj_test,paste0(outdir, 'rds/',id,'_ProjPC2Adj_TestSet.rds'))

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


