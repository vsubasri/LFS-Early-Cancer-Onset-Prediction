library(ggplot2)
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(ggpubr)

# Working directory should be project root when called from shell scripts
data <- readRDS("data/rds/NoobCorrected_beta_ProjPC2Adj_lfs_3UTR.rds")

#################################################
######### Plot wt tumor type breakdown ##########
#################################################

wt_clin <- read.csv("scripts/pipeline/Resources/wt_samples.csv")
ggplot(wt_clin,aes(x=cancer_diagnosis, y=..count..)) +
  geom_bar(position="dodge",color="black",fill="grey",width=0.85) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y= "Number of patients", x = "Cancer type") 

#######################################
######### Plot p53 Breakdown ##########
#######################################

p53<- data.frame(
  group=c("Mut (n=497)","Wt (n=132)"),
  value=c(497,132)
)

ggplot(p53, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(label = group), position = position_stack(vjust = 0.5),color = "black", size=4) +
  scale_fill_brewer(palette="Set1")


sex<- data.frame(
  group=c("Male","Female"),
  value=c(157,222)
)

#######################################
######### Plot Sex Breakdown ##########
#######################################
ggplot(sex, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(label = group), position = position_stack(vjust = 0.5),color = "black", size=6) +
  scale_fill_brewer(palette="Set2")

##########################################
######### Plot Family Breakdown ##########
##########################################

family <- as.data.frame(table(data$family)) %>% dplyr::group_by(Freq) %>% dplyr::count()
family$Freq <- as.factor(family$Freq)
ggplot(family,aes(x=Freq, y=n)) +
  geom_bar(stat="identity",color="black",fill="grey") + 
  theme_minimal() + 
  labs(y= "Number of Families", x = "# of Individuals within each Family \n with Methylation Profiling") 

###########################################
######### Plot Features by Chrom ##########
###########################################

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% data.frame()
feats <- read.csv(feature_path, sep="\t", row.names=NULL)
feats$chr <- factor(anno$chr[match(feats$probe,anno$Name)],levels=c(paste("chr",1:22,sep=""),"chrX","chrY","chrM"))

ggplot(feats,aes(x=chr, y=..count..)) +
  geom_bar(position="dodge",color="black",fill="grey") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y= "Number of Probes in Model", x = "Chromosome") 

##########################################################
######### Plot Age of Onset by LEF1 methylation ##########
##########################################################

data$earlyonset <- ifelse(is.na(data$ageofonset),"Cancer-free",ifelse(data$ageofonset < 72, "Early Onset","Late Onset"))
ggplot(data[c("LEF1","ageofonset","agesamplecollection")], aes(x=ageofonset, y=LEF1,color=agesamplecollection)) +
  geom_point() +
  stat_cor(method="pearson") + 
  theme_minimal()

##########################################################
######### Plot Age of Onset by RET methylation ###########
##########################################################

ggplot(data[c("RET","earlyonset")], aes(x=earlyonset, y=TSR1)) +
  geom_boxplot()  + 
  stat_compare_means(comparisons = list(c("Cancer-free", "Early Onset"),c("Late Onset","Early Onset") ), label="p.format",method="wilcox.test")


##########################################################

cell_props <- readRDS('~/Documents/lfs_ageofonset/cell_props/idol_cell_props_lfs.rds')
cell_props$sample_timing <- ifelse(cell_props$cancer1_age_diff > 0 | is.na(cell_props$cancer1_age_diff),"After","Before")
data <- cell_props[!duplicated(cell_props$ids),]
data <- data[!duplicated(data$tm_donor),]
data <- data[!is.na(data$agesamplecollection),]
data <- data[!is.na(data$ageofonset),]
data$early_onset <- ifelse(data$ageofonset < 72, "Early","After")
counts <- ddply(data, .(data$tissue_type, data$sample_timing), nrow)
colnames(counts) <- c("tissue_type","sample_timing","count")
counts <- counts[counts$tissue_type != "Unaffected",]
ggplot(counts,aes(tissue_type,count,fill=sample_timing)) +
  geom_bar(stat = "identity",position="dodge", colour="black",alpha=0.8) + 
  theme_minimal() + 
  xlab("Cancer Types") +
  ylab("Number of LFS Patients")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

tissue_ratio <- reshape(counts, idvar = "tissue_type", timevar = "sample_timing", direction = "wide")
tissue_ratio$ratio <- tissue_ratio$count.Before/tissue_ratio$count.After
tissue_ratio$`Sample Count` <- tissue_ratio$count.After+tissue_ratio$count.Before

tissue_plot <- ggplot(tissue_ratio,aes(tissue_type,ratio,fill=`Sample Count`)) +
  geom_bar(stat = "identity",position="dodge", colour="black",alpha=0.8) +
  theme_minimal() +
  xlab("Cancer Type") +
  ylab("Before/After")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank(),axis.title.y=element_blank()) +
  ylim(0, 0.7)

age_counts <- ddply(data, .(data$early_onset, data$sample_after), nrow)
colnames(age_counts) <- c("early_onset","sample_timing","count")
age_ratio <- reshape(age_counts, idvar = "early_onset", timevar = "sample_timing", direction = "wide")
age_ratio$ratio <- age_ratio$count.Before/age_ratio$count.After
age_ratio$`Sample Count` <- age_ratio$count.After+age_ratio$count.Before

age_plot <- ggplot(age_ratio,aes(early_onset,ratio,fill=`Sample Count`)) +
  geom_bar(stat = "identity",position="dodge", colour="black",alpha=0.8) + 
  theme_minimal() + 
  xlab("Cancer Onset") +
  ylab("Before/After")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylim(0, 0.7)

ggarrange(age_plot, tissue_plot, align = "hv",widths = c(0.3, 1),common.legend = TRUE,legend="bottom")


