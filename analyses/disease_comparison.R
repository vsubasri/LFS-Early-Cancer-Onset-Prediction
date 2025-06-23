library(ggplot2)
library(stringr)
library(ggpubr)
library(umap)

# Working directory should be project root when called from shell scripts

feats <- read.csv("data/features.txt", sep="\t", row.names=NULL)
all <- readRDS(file.path(dir,'data/beta_all.rds'))
all <- cbind(do.call(rbind,str_split(rownames(all),'_'))[,-1],all)
ms <- read.csv('data/GSE_metadata/GSE43976_metadata.txt',sep='\t')
ra <- read.csv('data/GSE_metadata/GSE42861_metadata.txt',sep='\t')
ms2 <- read.csv('data/GSE_metadata/GSE106648_metadata.txt',sep='\t')
ms3 <- read.csv('data/GSE_metadata/GSE130030_metadata.txt',sep='\t')
breast <- read.csv('data/GSE_metadata/GSE51057_metadata.txt',sep='\t')
colon <- read.csv('data/GSE_metadata/GSE51032_metadata.txt',sep='\t')
all <- all[!all$`1` %in% ms$X.Sample_geo_accession[ms$X.Sample_source_name_ch1 != "Peripheral blood"],]
all <- all[!all$`1` %in% ms3$X.Sample_geo_accession,]
all <- all[all$array == "450",]
lfs <- readRDS("data/rds/NoobCorrected_beta_ProjPC2Adj_lfs_3UTR.rds")
all$`1` <- str_replace(all$`1`,"X","")
all$sentrixid <- paste0(all$`1`,'_',all$`2`)
all$beta_id <- ifelse(all$sentrixid %in% as.character(na.omit(unique(lfs$SentrixID[lfs$ageofonset < 72]))),"lfsMut (Early Onset)",all$beta_id)
mut_exclude_ids <- lfs$SentrixID[(is.na(lfs$ageofonset) & lfs$cancer_diagnosis != "Unaffected")]
all <- all[!all$sentrixid %in% mut_exclude_ids,]
all$beta_id <- ifelse(all$beta_id == "lfsMut","lfsMut (Late Onset/Cancer-free)",all$beta_id)
wt_lfs <- read.csv('wt_samples.csv')
exclude_ids <- c("5760666053_R02C02","5760666057_R02C01","5760666057_R06C02", "5771955054_R06C01", "5760666003_R06C02","5760666035_R06C02","5760666022_R05C01","5771955054_R06C02")
wt_exclude_ids <- wt_lfs$SentrixID[(is.na(wt_lfs$age_diagnosis) & wt_lfs$cancer_diagnosis != "Unaffected")]
exclude_ids <- c(wt_exclude_ids,exclude_ids)
all <- all[!all$sentrixid %in% exclude_ids,]
all$beta_id <- ifelse(all$sentrixid %in% as.character(na.omit(unique(wt_lfs$SentrixID[wt_lfs$age_diagnosis < 72]))),"lfsWt (Early Onset)",all$beta_id)
all$beta_id <- ifelse(all$beta_id == "lfsWt","lfsWt (Late Onset/Cancer-free)",all$beta_id)
#all$beta_id <- ifelse(all$sentrixid %in% wt_lfs$SentrixID[wt_lfs$cancer_diagnosis_diagnoses == "Unaffected"],"lfsWt (Cancer-Free)",all$beta_id)
all$beta_id <- str_replace(all$beta_id,"breast","Breast")
all$beta_id <- str_replace(all$beta_id,"colon","Colon")
all$beta_id <- str_replace(all$beta_id,"control","Control")
all$beta_id <- str_replace(all$beta_id,"ms","Multiple Sclerosis")
all$beta_id <- str_replace(all$beta_id,"ra","Rheumatoid Arthritis")

final_meth <- all[6:(length(all)-2)]
u <- umap(final_meth,n_neighbors=50)
u <- cbind(data.frame(u$layout))
u$dataset <- all$beta_id
u$array <- all$array

ggplot(u, aes(X1,X2,color=dataset)) +
  geom_point(alpha=0.4) + 
  theme_minimal() + 
  scale_colour_manual(values = c( "green","orange","light blue", "red", "dark blue", "purple", "pink","brown","black")) +
  guides(color=guide_legend(title="PBL Methylation Dataset"))

ggplot(u[u$dataset %in% c("lfsMut (Late Onset/Cancer-free)","lfsMut (Early Onset)","lfsWt (Late Onset/Cancer-free)","lfsWt (Early Onset)"),], aes(X1,X2,color=dataset)) +
  geom_point() + 
  theme_minimal() + 
  scale_colour_manual(values = c( "red", "dark blue", "purple","pink")) +
  guides(color=guide_legend(title="PBL Methylation Dataset"))

u$tp53_status <- ifelse(grepl("lfsMut",u$dataset),"LFS (Mut)", ifelse(grepl("lfsWt",u$dataset),"LFS (Wt)", u$dataset))
u$onset_status <- ifelse(grepl("Early",u$dataset),"Early Onset", ifelse(grepl("Late",u$dataset),"Late Onset/Cancer-Free", u$dataset))
ggplot(u[u$dataset %in% c("lfsMut (Late Onset/Cancer-free)","lfsMut (Early Onset)","lfsWt (Late Onset/Cancer-free)","lfsWt (Early Onset)"),], aes(onset_status,X1,fill=onset_status)) +
  geom_boxplot() + 
  theme_minimal() + 
  facet_wrap(~tp53_status)+
  stat_compare_means(comparisons = list(c("Early Onset", "Late Onset/Cancer-Free"))) + 
  theme(legend.position = "bottom")

ggplot(u[u$dataset %in% c("lfsMut (Late Onset/Cancer-free)","lfsMut (Early Onset)","lfsWt (Late Onset/Cancer-free)","lfsWt (Early Onset)"),], aes(dataset,X2,fill=dataset)) +
  geom_boxplot() + 
  theme_minimal() + 
  stat_compare_means(comparisons = list(c("lfsMut (Early Onset)", "lfsMut (Late Onset/Cancer-free)"),
                                        c("lfsWt (Early Onset)", "lfsWt (Late Onset/Cancer-free)")))


p <- prcomp(final_meth,scale=TRUE)
p <- cbind(data.frame(p$x))
p$dataset <- all$beta_id
p$ids <- all$ids
p$array <- all$array

ggplot(p, aes(PC1,PC2,color=dataset)) +
  geom_point(alpha = 0.2) + 
  theme_minimal()  

