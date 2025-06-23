## load up the packages we will need:

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(umap))
suppressMessages(library(argparse))

## set up parser

parser <- ArgumentParser()
parser$add_argument("--id", action="store")
parser$add_argument("--seed", action="store")
parser$add_argument("--nsplit", action="store",type="integer")
parser$add_argument("--outdir", action="store")

args <- parser$parse_args()
id <- args$id
seed <- args$seed
nsplit <- args$nsplit
outdir <- args$outdir

# Working directory should be project root when called from shell scripts
source('preprocessing/util_functions.R')

set.seed(seed)

cat("[ Generate train/test/val split  ]","\n")
data <- readRDS(paste0(outdir,"rds/",id,".rds"))
data <- data[!is.na(data$ids),]
data <- data[!data$Meth %in% c("Problem"),]

ind <- grep("cg", colnames(data))[1]-1

datasplits <- getTrainTestVal(data,seed,nsplit)
data_train <- datasplits[[1]]
data_test <- datasplits[[2]]
data_val <- datasplits[[3]]
rm(data)

p_train <- prcomp(data_train[(ind+1):length(data_train)],scale=T)
saveRDS(p_train,paste0('rds/',id,'_',seed,'_batchEffectRemoval.rds'))
p_train_clin <- cbind(data_train[1:ind],p_train$x)

cat("[ Find PC most correlated to plate ]",'\n')
# Find pcs correlated with age of sample collection
pval <- list()
for (i in 1:dim(p_train$x)[1]) {
  aov_pc <- aov(p_train$x[, i] ~ p_train_clin$Project)
  pval[[i]] <- summary(aov(p_train$x[, i] ~ p_train_clin$Project))[[1]][[1,"Pr(>F)"]]
}

pca_corr <- data.frame(pval=do.call('rbind', pval))
pca_corr$rank <- rank(pca_corr$pval)
pca_corr$pc <- paste0("PC",row.names(pca_corr))

maxpc <- pca_corr$pc[pca_corr$rank == 1] ; maxpc_p <- pca_corr$pval[pca_corr$rank == 1]
maxpc2 <- pca_corr$pc[pca_corr$rank == 2] ; maxpc_p2 <- pca_corr$pval[pca_corr$rank == 2]
cat(paste0("[ Max batch-associated PC ] :",maxpc,"\t","[ Max PC Association with Batch ] :",maxpc_p,"\n"))
cat(paste0("[ Max batch-associated PC 2 ] :",maxpc2,"\t","[ Max PC 2 Association with Batch ] :",maxpc_p2,"\n"))

cat("[ Remove PC most correlated with Batch Effect ]",'\n')
#https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com
Xhat <- p_train$x[, !(colnames(p_train$x) %in% c(maxpc))] %*% t(p_train$rotation[, !(colnames(p_train$rotation) %in% c(maxpc))])
beta_adj <- scale(Xhat, center = -(colMeans(data_train[(ind+1):length(data_train)])), scale = T)
data_adj <- cbind(data_train[1:ind],beta_adj)
saveRDS(data_adj,paste0(outdir, 'rds/',id,'_ProjPC2Adj_TrainingSet.rds'))

## Plot training set correction ##
cat("[ Plotting PCA of corrected trainig data ]","\n")
pc <- prcomp(beta_adj,scale=T)
pc_clin <- cbind(data_adj[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir, 'Output/',id,'_ProjPC2Adj_TrainingSet_PCA.csv'),quote=F)

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

