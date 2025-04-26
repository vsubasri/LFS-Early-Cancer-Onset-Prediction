library(ggplot2)
library(ggpubr)

train_dat <- readRDS('~/research/lfs_methylation/Data_objects/NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_TSS200_TrainingSet.rds') 
test_dat <- readRDS('~/research/lfs_methylation/Data_objects/NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_TSS200_TestSet.rds') 
val_dat <- readRDS('~/research/lfs_methylation/Data_objects/NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_TSS200_ValidationSet.rds') 

remove_outliers <- function(pc,n) {
  pc$PC1 <- scale(pc$PC1)
  pc$PC2 <- scale(pc$PC2)
  u_thres_pc1 <- mean(pc$PC1) + n*sd(pc$PC1)
  l_thres_pc1 <- mean(pc$PC1) - n*sd(pc$PC1)
  u_thres_pc2 <- mean(pc$PC2) + n*sd(pc$PC2)
  l_thres_pc2 <- mean(pc$PC2) - n*sd(pc$PC2)
  pc2 <- pc[pc$PC1 < u_thres_pc1 & pc$PC1 > l_thres_pc1, ]
  pc2 <- pc2[pc2$PC2 < u_thres_pc2 & pc2$PC2 > l_thres_pc2, ]
  outlier_id <- as.character(pc$SentrixID[!pc$SentrixID %in% pc2$SentrixID])
  cat(paste0(length(outlier_id)," outliers removed","\n"))
  keep_id <- as.character(pc2$SentrixID)
  return(keep_id)
}

pc <- read.csv('~/research/lfs_methylation/Output/Seed9876_S40/Noob_beta_PCA.csv')
pc$Project <- factor(pc$Project,levels=c("SUB12648", "SUB14749", "SUB14750", "SUB14751", "TORONTO", "MONTREAL"))
pc$array <- factor(pc$array)
pc_corrected <- read.csv('~/research/lfs_methylation/Output/Seed9876_S40/NoobCorrected_after_covs_beta_ProjPC2Adj_PCA.csv')
pc_corrected$Project <- factor(pc_corrected$Project,levels=c("SUB12648", "SUB14749", "SUB14750", "SUB14751", "TORONTO", "MONTREAL"))
pc_corrected$array <- factor(pc_corrected$array)

library(RColorBrewer)

plot_pca <- function(data)  {
  data <- data[data$SentrixID %in% remove_outliers(data,3),]
  cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$Project)))
  plot <- ggplot(data, aes(PC1, PC2, color = Project,shape=array)) +
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal(base_size = 18)
  return(plot)
}

train_plot <- plot_pca(pc[pc$ids %in% train_dat$ids,])
test_plot <- plot_pca(pc[pc$ids %in% test_dat$ids,])
val_plot <- plot_pca(pc[pc$ids %in% val_dat$ids,])

ctrain_plot <- plot_pca(pc_corrected[pc_corrected$ids %in% train_dat$ids,])
ctest_plot <- plot_pca(pc_corrected[pc_corrected$ids %in% test_dat$ids,])
cval_plot <- plot_pca(pc_corrected[pc_corrected$ids %in% val_dat$ids,])

pca_plots <- ggarrange(train_plot, test_plot, val_plot,
                       ctrain_plot,ctest_plot,cval_plot,
                                 labels = c("Train (before)", "Validation (before)", "Test (before)",
                                            "Train (after)", "Validation (after)", "Test (after)"),
                                 ncol = 3,nrow=2,common.legend = TRUE)
pca_plots

