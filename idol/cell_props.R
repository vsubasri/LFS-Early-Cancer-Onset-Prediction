library(ggplot2)
library(umap)
library(corrplot)
library(umap)
library(ggsci)
library(reshape2)
library(ggpubr)

# Working directory should be project root when called from shell scripts
data <- readRDS("data/rds/NoobCorrected_beta_ProjPC2Adj_lfs_3UTR.rds")
bestmodel <- "checkpoint/NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_tp53_72_svmRadial_ageofonset_results"

# Load checkpoint model results  
checkpoint_data <- readRDS(paste0(bestmodel, '.rds'))
ROCobj_val <- checkpoint_data$ROCInfoVal
val_results <- ROCobj_val[[1]]
ROCobj_test <- checkpoint_data$ROCInfoTest 
test_results <- ROCobj_test[[1]]

data$Dataset <- ifelse(data$ids %in% val_results$ids, "Validation", ifelse(data$ids %in% test_results$ids,"Test","Training"))

cell_props <- readRDS('~/Documents/lfs_ageofonset/cell_props/idol_cell_props_lfs.rds')
cell_props$sample_timing <- ifelse(cell_props$cancer1_age_diff <= 0,"Before","After")
cell_props$sample_timing <- ifelse(cell_props$cancer_diagnosis == "Unaffected","Cancer-Free",cell_props$sample_timing)

u <- umap(cell_props[c("CD8T","CD4T","NK","Bcell","Mono","Neu")])
u <- cbind(cell_props,data.frame(u$layout))
u$earlyonset <- ifelse((u$ageofonset > 72 | is.na(u$ageofonset)),"N","P")
cutoff <- 0.4
val_FN <- val_results$ids[val_results$test_label == 1 & val_results$test_pred_calibrated.Yes < cutoff]
test_FN <- test_results$ids[test_results$test_label == 1 & test_results$test_pred_calibrated.Yes < cutoff]
val_FP <- val_results$ids[val_results$test_label == 0 & val_results$test_pred_calibrated.Yes > cutoff]
test_FP <- test_results$ids[test_results$test_label == 0 & test_results$test_pred_calibrated.Yes > cutoff]
u$false <- ifelse(u$ids %in% c(val_FN,test_FN,val_FP,test_FP),"F","T")

results <- rbind(val_results, test_results)
pred_u <- u[u$ids %in% results$ids,]
pred_u$prob <- results$test_pred_calibrated.Yes[match(pred_u$ids, results$ids)]
res <- round(cor(pred_u[c("CD8T","CD4T","NK","Bcell","Mono","Neu","prob")]),2)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) 

p <- prcomp(cell_props[c("CD8T","CD4T","NK","Bcell","Mono","Neu")], scale=TRUE)
summary(p)
p <- cbind(data.frame(p$x),cell_props)
p$earlyonset <- ifelse((p$ageofonset > 72 | is.na(p$ageofonset)),"N","P")
p$false <- ifelse(p$ids %in% c(val_FN,test_FN,val_FP,test_FP),"F","T")

ggplot(p, aes(PC1,PC2,color=interaction(earlyonset,false))) +
  geom_point() + 
  theme_minimal()

ggplot(pred_u, aes(X1,X2,color=interaction(earlyonset,false))) +
  geom_point(size=2,alpha=0.8) + 
  theme_minimal()

ggplot(pred_u, aes(Bcell, prob,color=interaction(false,earlyonset))) +
  geom_point(size=3) + 
  theme_classic() + 
  theme(axis.ticks.x=element_blank()) + 

plot_cellprops <- melt(pred_u[c("CD8T","CD4T","NK","Bcell","Mono","Neu","earlyonset","agesamplecollection","sample_timing","false")],id=c("earlyonset","sample_timing","agesamplecollection"),value=c("CD8T","CD4T","NK","Bcell","Mono","Neu"))
plot_cellprops$value <- as.numeric(plot_cellprops$value)
plot_cellprops$earlyonset <- ifelse(plot_cellprops$earlyonset=="N","Cancer after 6/Cancer-free","Cancer before 6")

data_adjusted <- compositions::acomp(pred_u[2:7] + 1e-6)
clr_transformed <- compositions::clr(data_adjusted)

# Convert clr_transformed to a data frame and prepare for ggplot
clr_df <- as.data.frame(clr_transformed)
clr_df$earlyonset <- pred_u$earlyonset
clr_df$sample_timing <- pred_u$sample_timing
clr_df$agesamplecollection <- pred_u$agesamplecollection

# Melting the data
plot_cellprops <- melt(clr_df, id.vars = c("earlyonset", "sample_timing", "agesamplecollection"))

# Convert factors if necessary
plot_cellprops$earlyonset <- ifelse(plot_cellprops$earlyonset == "N", "Cancer after 6/Cancer-free", "Cancer before 6")
plot_cellprops$value <- as.numeric(plot_cellprops$value)

# Plotting
ggplot(plot_cellprops, aes(x = variable, y = value, fill = earlyonset)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.5, size = 1.5) +
  scale_fill_brewer(palette = "Set1", name = "Group") +
  labs(
    x = "Cell Type",
    y = "Center Log Ratio (CLR) Transformed Proportion",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  stat_compare_means(aes(label = paste0("p=", ..p.adj..)),
                     method = "wilcox.test",
                     label.y = max(plot_cellprops$value) * 1.1,
                     size = 3.5,
                     tip.length = 0) +  # Remove connector line for cleaner look
  scale_y_continuous(labels = scales::comma)  # Format y-axis labels for better readability
#  scale_color_distiller(palette = 'Spectral')
#  geom_point(aes(color = sample_timing, group = earlyonset), position = position_dodge(width = 0.75), alpha = 0.5) +
 # scale_color_manual(values = c("Before" = "red", "After" = "black","Cancer-Free"="blue")) +

pred_u_adjusted <- compositions::acomp(pred_u[c("CD8T","CD4T","NK","Bcell","Mono","Neu")] + 1e-6)
pred_u_clr_transformed <- as.data.frame(compositions::clr(pred_u_adjusted))
pred_u_clr_transformed$prob <- pred_u$prob
res <- round(cor(pred_u_clr_transformed),2)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) 

library(RColorBrewer)
n <- length(names(table(u$cancer_diagnosis))[table(u$cancer_diagnosis)>4])
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(u[u$cancer_diagnosis %in% names(table(u$cancer_diagnosis))[table(u$cancer_diagnosis)>4],], aes(X1,X2,color=cancer_diagnosis)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_manual(values=col_vector)

n <- length(names(table(u$tissue_type))[table(u$tissue_type)>1])
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

tissue_types=names(table(u$tissue_type))[table(u$tissue_type)>1]
ggplot(u[u$tissue_type %in% tissue_types,], aes(X1,X2,color=tissue_type)) +
  geom_point(size=2) + 
  theme_minimal() + 
  scale_color_manual(values=col_vector)

library(scales)

ggplot(u, aes(X1,X2,color=ageofonset)) +
  geom_point(size=2) + 
  theme_minimal() + 
  scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))

ggplot(u, aes(X1,X2,color=agesamplecollection)) +
  geom_point(size=2) + 
  theme_minimal() + 
  scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))

ggplot(u, aes(X1,X2,color=interaction(false,earlyonset))) +
  geom_point(size=2,alpha=0.8) + 
  theme_minimal()

ggplot(u, aes(X1,X2,color=CD4T,shape=interaction(false))) +
  geom_point(size=2) + 
  theme_minimal() + 
  scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))

ggplot(u, aes(X1,X2,color=CD8T)) +
  geom_point(size=2) + 
  theme_minimal() + 
  scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))

ggplot(u, aes(X1,X2,color=Neu)) +
  geom_point(size=2) + 
  theme_minimal() + 
  scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))

ggplot(u, aes(X1,X2,color=Bcell)) +
  geom_point(size=2) + 
  theme_minimal() + 
  scale_colour_gradientn(colours = c("darkred", "orange", "yellow", "white"))

cellprops_all <- readRDS('~/Documents/lfs_ageofonset/prop_all_formatted.rds')
cellprops_all <- cellprops_all[cellprops_all$id %in% c("LFSMut","LFSWt"),]
cellprops_all <- melt(cellprops_all,value=c("CD8T","CD4T","NK","Bcell","Mono","Neu"))

ggplot(cellprops_all[cellprops_all$array == "450",], aes(variable, value,fill=id)) +
  geom_point(position=position_jitterdodge(),alpha=0.5) +
  geom_boxplot(outlier.shape=NA) + 
  theme_classic() + 
  theme(axis.ticks.x=element_blank()) + 
  xlab("Cell Type") + ylab("Proportion") + 
  stat_compare_means(aes(label=paste0("p=",..p.adj..)),label.y = 1.1, size = 3, method = "wilcox.test") +
  theme_classic() +
  scale_fill_jco() 
