library(dplyr)
library(ggplot2)

# Working directory should be project root when called from shell scripts

# load immune cell counts
#c <- readRDS("count_450.rds") # IDOL with 450k array
c <- readRDS("count_EPIC.rds") # IDOL wiht EPIC array (recommended)
counts <- c$counts
prop <- prop.table(counts, margin=1)
data <- as.data.frame(prop) # dataframe of proportion of each immune cell per sample

# load clincial data
clin <- read.csv("lfs_mut_clinical_comprehensive.csv")

# load sample sheets to correspond count data to clinical
schiffman <- read.csv("SampleSheet_schiffman.csv")
batch5 <- read.csv("SampleSheet_batch5.csv")
montreal <- read.csv("SampleSheet_montreal.csv")
toronto <- read.csv("SampleSheet_toronto.csv")

# dataframes with methylation ids and sample names
schiffman_id<- data.frame(SentrixID=paste0(schiffman$Sentrix_ID,"_",schiffman$Sentrix_Position),
                          ids=cbind(lapply(schiffman$Sample_Name, function(x) x[length(x)])))
batch5_id<- data.frame(SentrixID=paste0(batch5$Sentrix_ID,"_",batch5$Sentrix_Position),
                       ids=cbind(lapply(batch5$Sample_Name, function(x) x[length(x)])))
montreal_id<- data.frame(SentrixID=paste0(montreal$Sentrix_ID,"_",montreal$Sentrix_Position),
                         ids=cbind(lapply(montreal$Sample_Name, function(x) x[length(x)])))
toronto_id<- data.frame(SentrixID=paste0(toronto$Sentrix_ID,"_",toronto$Sentrix_Position),
                        ids=cbind(lapply(strsplit(toronto$Sample_Name,'#'), function(x) x[length(x)])))
methid2samplename <- bind_rows(schiffman_id, batch5_id, montreal_id, toronto_id)
methid2samplename$ids2 <- unlist(methid2samplename$ids)
methid2samplename$ids <- NULL

# load sample meth info
samples_450 <- readRDS("450.rds")
samples_850 <- readRDS("850.rds")
samples <- rbind(samples_450, samples_850)
samples$data_sample <- basename(samples$Basename) # sample name corresponding to row names of data dataframe

# merge data together
ids_clin <- merge(methid2samplename, clin, by.x="ids2", by.y="sample", all.x=TRUE)
ids_clin_samples <- merge(samples, ids_clin, by.x="data_sample", by.y="SentrixID", all.x=TRUE)
all <- merge(ids_clin_samples, data, by.x="data_sample", by.y="row.names", all.x=TRUE)


### PLOT CONTINUOUS ###
# plot two continuous variables
plot_reg <- function(data, type, x, y){
  v1 <- data[[x]]
  v2 <- data[[y]]
  m_x <- max(v1, na.rm=TRUE) # max used for placement of text
  m_y <- max(v2, na.rm=TRUE) # max used for placement of text
  model <- lm(v2~v1, data=data) # linear model
  r_sqr <- round(summary(model)$r.square, 2) # r square
  change <- round(summary(model)$coefficients[2,1], 4) # estimate from lm 
  p <- round(summary(model)$coefficients[2,4], 4) # p value from lm
  annot <- paste("R^2 =", r_sqr, "\n", "Estimate =", change, "\n", "p value =", p) # annotation string
  file <- paste0("figures/", type, "_", x, "_", y, ".pdf") 
  pdf(file)
  
  g <- ggplot(data) +
    aes(x = v1, y = v2) +
    geom_point(colour = "#0c4c8a") +
    geom_smooth(method="lm", colour = "red") + # line on plot
    theme_bw() +
    xlab(x) + 
    ylab(y) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    annotate(x=m_x, y=m_y-100, # annotation
             label= annot, 
             geom="text", size=5, hjust=1)
  show(g)
  dev.off()
}

# age of onset
# estimate is 1 change in x results in that change in y
# p value null hypothesis is coefficient is 0 (no relationship b/w x and y)
plot_reg(all, "850", "CD8T", "ageofonset")
plot_reg(all, "850","CD4T", "ageofonset")
plot_reg(all, "850","NK", "ageofonset")
plot_reg(all, "850","Bcell", "ageofonset")
plot_reg(all, "850","Mono", "ageofonset")
plot_reg(all, "850","Neu", "ageofonset")

# age of onset
plot_reg(all, "850","CD8T", "agesamplecollection")
plot_reg(all, "850","CD4T", "agesamplecollection")
plot_reg(all, "850","NK", "agesamplecollection")
plot_reg(all, "850","Bcell", "agesamplecollection")
plot_reg(all, "850","Mono", "agesamplecollection")
plot_reg(all, "850","Neu", "agesamplecollection")

# multivariate linear regression
# age of onset
model <- lm(ageofonset~CD8T+CD8T+NK+Bcell+Mono+Neu, data=all)
summary(model)

# age sample collection
model <- lm(agesamplecollection~CD8T+CD8T+NK+Bcell+Mono+Neu, data=all)
summary(model)


### PLOT CATEGORICAL ###
# cancer status
plot_cat <- function(data, type, x, y){
  temp <- data %>% 
    drop_na({{x}}) %>% # remove rows with NA in column of interest
    group_by(across({{x}})) %>%
    filter(n()>4) %>% # must have 5 or more samples to graph
    as.data.frame()
  v1 <- temp[[x]]
  v2 <- temp[[y]]
  combo <- combn(unique(v1), 2) # combinations for pairwise statistical test
  comp <- as.list(as.data.frame(combo))
  filepath <- paste0("figures/", type, "_", x, "_", y, ".pdf") 
  #pdf(filepath)
  
  g <- ggplot(temp) + 
    aes_string(x={{x}}, y={{y}}, fill={{x}}) + 
    geom_boxplot(notch=FALSE) + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    coord_cartesian(clip = "off") + # allows text outside plot
    stat_compare_means() # wilcox if 2 or kruskal if 3+
  show(g)
  #dev.off()
}


# cancer at draw
all$cancer_atdraw_fac <- ifelse(all$cancer_atdraw=="No", 0, 1) # affected or unaffected column
plot_cat(all, "850","cancer_atdraw", "CD8T")
plot_cat(all, "850","cancer_atdraw", "CD4T")
plot_cat(all, "850","cancer_atdraw", "NK")
plot_cat(all, "850","cancer_atdraw", "Bcell")
plot_cat(all, "850","cancer_atdraw", "Mono")
plot_cat(all, "850","cancer_atdraw", "Neu")

model <- glm(cancer_atdraw_fac~CD8T+CD8T+NK+Bcell+Mono+Neu, data=all)
summary(model)


# first cancer diagnosis
plot_cat(all, "850","cancer_diagnosis", "CD8T")
plot_cat(all, "850","cancer_diagnosis", "CD4T")
plot_cat(all, "850","cancer_diagnosis", "NK")
plot_cat(all, "850","cancer_diagnosis", "Bcell")
plot_cat(all, "850","cancer_diagnosis", "Mono")
plot_cat(all, "850","cancer_diagnosis", "Neu")


# affected or unaffected
all$affected <- ifelse(all$cancer_diagnosis=="Unaffected", 0, 1) # affected or unaffected column
plot_cat(all, "450","affected", "CD8T")
plot_cat(all, "450","affected", "CD4T")
plot_cat(all, "450","affected", "NK")
plot_cat(all, "450","affected", "Bcell")
plot_cat(all, "450","affected", "Mono")
plot_cat(all, "450","affected", "Neu")

# age of onset
model <- glm(affected~CD8T+CD8T+NK+Bcell+Mono+Neu, data=all)
summary(model)
