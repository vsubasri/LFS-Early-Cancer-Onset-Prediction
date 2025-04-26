library(ggplot2)
library(ggpubr)

data <- readRDS('/Users/vallijahsubasri/Documents/lfs_ageofonset/data/rds/NoobCorrected_beta_ProjPC2Adj.rds')
data <- data[!is.na(data$agesamplecollection), ]
age_probes <- read.csv('/Users/vallijahsubasri/Documents/lfs_ageofonset/scripts/pipeline/Resources/Horvath_PedBE_probes.csv',stringsAsFactors = F)
data_age <- data[age_probes$probe[age_probes$probe %in% colnames(data)]]
p_age <- prcomp(data_age, scale=TRUE)
summary(p_age)
plot_p_age <- cbind(data.frame(p_age$x),data[1:67])
plot_p_age$agesamplecollection_cutoff <- factor(ifelse(plot_p_age$agesamplecollection < 72,"Age of Sample Collection <= 6 Years", "Age of Sample Collection > 6 Years"))
c_pc1_age <- cor.test(plot_p_age$PC1,plot_p_age$agesamplecollection,method = "pearson")
age_pc1_label <- paste0("PC1: R = ",round(c_pc1_age$estimate,3), " , Variance Explained = ",round(summary(p_age)$importance[2,1]*100,1))
c_pc2_age <- cor.test(plot_p_age$PC2,plot_p_age$agesamplecollection)
age_pc2_label <- paste0("PC2: R = ",round(c_pc2_age$estimate,3), " , Variance Explained = ",round(summary(p_age)$importance[2,2]*100,2))
age_plot <- ggplot(plot_p_age, aes(PC1,PC2,color=agesamplecollection_cutoff)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_manual(values = c("Age of Sample Collection <= 6 Years" = "#2E86C1", 
                                "Age of Sample Collection > 6 Years" = "#E67E22")) +
  labs(color = "Cutoff") +
  annotate("text", size=3, x=0, y=-27, label=age_pc1_label) + 
  annotate("text", size=3, x=0, y=-30,label=age_pc2_label)
age_plot

age_corr <- list()
for (probe in colnames(data_age[68:length(data_age)])) {
  print(probe)
  age_corr[probe] <- cor.test(data_age[,probe],data$agesamplecollection)$estimate
}
age_corr <- data.frame(do.call('rbind',age_corr))
colnames(age_corr) <- c("corr")
age_corr$type <- "age-associated"

features <- read.csv('/Users/vallijahsubasri/Documents/lfs_ageofonset/scripts/pipeline/Resources/features.txt',sep='\t')
data_model <- data[features$probe]
p_model <- prcomp(data_model[68:length(data_model)], scale=TRUE)
summary(p_model)
plot_p_model  <- cbind(data.frame(p_model$x),data[1:67])
plot_p_model$agesamplecollection_cutoff <- factor(ifelse(plot_p_model$agesamplecollection < 72, "Age of Sample Collection <= 6 Years", "Age of Sample Collection > 6 Years"))
c_pc1_model <- cor.test(plot_p_model$PC1,plot_p_model$agesamplecollection)
model_pc1_label <- paste0("PC1: R = ",round(c_pc1_model$estimate,3), " , Variance Explained = ",round(summary(p_model)$importance[2,1]*100,1))
c_pc2_model <- cor.test(plot_p_model$PC2,plot_p_model$agesamplecollection)
model_pc2_label <- paste0("PC2: R = ",round(c_pc2_model$estimate,3), " , Variance Explained = ",round(summary(p_model)$importance[2,2]*100,1))

notage_plot <- ggplot(plot_p_model, aes(PC1,PC2,color=agesamplecollection_cutoff)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values = c("Age of Sample Collection <= 6 Years" = "#2E86C1", 
                                "Age of Sample Collection > 6 Years" = "#E67E22")) +
  labs(color = "Cutoff") +
  annotate("text", size=3, x=-8, y=-16, label=model_pc1_label) + 
  annotate("text", size=3, x=-8, y=-17.5,label=model_pc2_label)
notage_plot

model_corr <- list()
for (probe in colnames(data_model[68:length(data_model)])) {
  print(probe)
  model_corr[probe] <- cor.test(data_model[,probe],data$agesamplecollection)$estimate
}
model_corr <- data.frame(do.call('rbind',model_corr))
colnames(model_corr) <- c("corr")
model_corr$type <- "model"

ggarrange(age_plot,notage_plot,common.legend = TRUE,legend="bottom",
          labels=c("A","B"))

corr_dist <- rbind(model_corr,age_corr)
corr_dist$abs_corr <- abs(corr_dist$corr)

ggplot(corr_dist,aes(x=abs_corr,fill=type)) + 
  theme_minimal() +
  geom_density(alpha=0.5)


###############
data <- readRDS('/Users/vallijahsubasri/Documents/lfs_ageofonset/NoobCorrected_beta_ProjPC2Adj.rds')
age_probes <- read.csv('/Users/vallijahsubasri/Documents/lfs_ageofonset/Horvath_PedBE_probes.csv',stringsAsFactors = F)
data_notage <- data[colnames(data)[!colnames(data) %in% age_probes$probe]]
p_age <- prcomp(data, scale=TRUE)
plot_p_age <- cbind(data.frame(p_age$x),data[1:67])
c_pc1_age <- cor.test(plot_p_age$PC1,plot_p_age$agesamplecollection,method = "pearson")
age_pc1_label <- paste0("PC1: R = ",round(c_pc1_age$estimate,3), " , Variance Explained = ",round(summary(p_age)$importance[2,1]*100,1))
c_pc2_age <- cor.test(plot_p_age$PC2,plot_p_age$agesamplecollection)
age_pc2_label <- paste0("PC2: R = ",round(c_pc2_age$estimate,3), " , Variance Explained = ",round(summary(p_age)$importance[2,2]*100,2))


p_notage <- prcomp(data_notage, scale=TRUE)
plot_p_notage <- cbind(data.frame(p_notage$x),data[1:67])
c_pc1_notage <- cor.test(plot_p_notage$PC1,plot_p_notage$agesamplecollection,method = "pearson")
notage_pc1_label <- paste0("PC1: R = ",round(c_pc1_notage$estimate,3), " , Variance Explained = ",round(summary(p_notage)$importance[2,1]*100,1))
c_pc2_notage <- cor.test(plot_p_notage$PC2,plot_p_notage$agesamplecollection)
notage_pc2_label <- paste0("PC2: R = ",round(c_pc2_notage$estimate,3), " , Variance Explained = ",round(summary(p_notage)$importance[2,2]*100,2))


age_plot <- ggplot(plot_p_age, aes(PC1,PC2,color=agesamplecollection)) +
  geom_point() + 
  theme_minimal() + 
  annotate("text", size=3, x=0, y=-45, label=age_pc1_label) + 
  annotate("text", size=3, x=0, y=-49,label=age_pc2_label)
age_plot
