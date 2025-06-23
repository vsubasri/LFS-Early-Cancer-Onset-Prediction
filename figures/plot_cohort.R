library(ggplot2)
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)

# Working directory should be project root when called from shell scripts

bestmodel <- "checkpoint/NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_tp53_72_svmRadial_ageofonset_results"

# Load checkpoint model results  
checkpoint_data <- readRDS(paste0(bestmodel, '.rds'))
ROCobj_val <- checkpoint_data$ROCInfoVal
val_results <- ROCobj_val[[1]]
ROCobj_test <- checkpoint_data$ROCInfoTest 
test_results <- ROCobj_test[[1]]

data <- readRDS("data/rds/Noob_beta.rds")
data <- data[!duplicated(data$tm_donor),]
#data <- readRDS("data/rds/NoobCorrected_beta_ProjPC2Adj_lfs_3UTR.rds")
data <- data[!data$ids %in% c(2765,2957,3365,3367,3503,3885,4033,3301,3692,3733,3884,3886,3887,3888,3889,3892,4036,4037,4164,4165,4169,3304),]
data$Dataset <- ifelse(data$ids %in% val_results$ids, "Validation", ifelse(data$ids %in% test_results$ids,"Test","Training"))

wt_clin <- read.csv("scripts/pipeline/Resources/wt_samples.csv")
wt_clin$tm_donor <- as.character(wt_clin$tm_donor)
wt_clin$age_diagnosis <- wt_clin$age_diagnosis/12
wt_clin$age_sample_collection <- wt_clin$age_diagnosis/12

# Prepare WT data
wt_data <- wt_clin %>%
  select(ids, tm_donor, age_diagnosis, age_sample_collection, cancer_diagnosis, Project, SentrixID,p53_germline) %>%
  mutate(
    ids=ids,
    tm_donor = tm_donor,
    Project = Project, 
    SentrixID = SentrixID,
    Dataset = "WT Control",
    ageofonset = age_diagnosis,
    agesamplecollection = age_sample_collection,
    cancer_diagnosis = cancer_diagnosis,
    p53=p53_germline
  ) %>%
  select(ids, tm_donor, ageofonset, agesamplecollection, cancer_diagnosis, Project, Dataset, SentrixID,p53) 

# Combine LFS and WT data
combined_data <- data %>%
  select(ids, tm_donor, Project, SentrixID, p53, Dataset, ageofonset, agesamplecollection, cancer_diagnosis) %>%
  bind_rows(wt_data)

write.csv(combined_data,'Supplementary Table 1.csv',row.names = F)

# Perform Wilcoxon rank-sum test
wilcox_result <- wilcox.test(agesamplecollection ~ cancerstatus, data = data)

# Print results
print("Age of Sample Collection: Cancer vs Unaffected")
print(sprintf("Wilcoxon rank-sum test p-value: %.2e", wilcox_result$p.value))

# Get summary statistics - alternative version without n()
# Calculate confidence intervals and range for each group
summary_stats <- data.frame(
  cancerstatus = c("Cancer", "Unaffected"),
  count = c(sum(data$cancerstatus == "Cancer"), 
            sum(data$cancerstatus == "Unaffected")),
  median_age = c(median(data$agesamplecollection[data$cancerstatus == "Cancer"]/12, na.rm = TRUE),
                 median(data$agesamplecollection[data$cancerstatus == "Unaffected"]/12, na.rm = TRUE)),
  mean_age = c(mean(data$agesamplecollection[data$cancerstatus == "Cancer"]/12, na.rm = TRUE),
               mean(data$agesamplecollection[data$cancerstatus == "Unaffected"]/12, na.rm = TRUE)),
  sd_age = c(sd(data$agesamplecollection[data$cancerstatus == "Cancer"]/12, na.rm = TRUE),
             sd(data$agesamplecollection[data$cancerstatus == "Unaffected"]/12, na.rm = TRUE)),
  min_age = c(min(data$agesamplecollection[data$cancerstatus == "Cancer"]/12, na.rm = TRUE),
              min(data$agesamplecollection[data$cancerstatus == "Unaffected"]/12, na.rm = TRUE)),
  max_age = c(max(data$agesamplecollection[data$cancerstatus == "Cancer"]/12, na.rm = TRUE),
              max(data$agesamplecollection[data$cancerstatus == "Unaffected"]/12, na.rm = TRUE))
)

# Add confidence intervals
summary_stats$ci_lower <- with(summary_stats, 
                              mean_age - qt(0.975, count-1) * (sd_age/sqrt(count)))
summary_stats$ci_upper <- with(summary_stats, 
                              mean_age + qt(0.975, count-1) * (sd_age/sqrt(count)))

print(summary_stats)

# Add boxplot to visualize the difference
ggplot(data, aes(x = cancerstatus, y = agesamplecollection/12)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Cancer Status", 
       y = "Age of Sample Collection (years)",
       title = sprintf("p = %.2e", wilcox_result$p.value))


#######################################
######### Figure 2B (cohort) ##########
#######################################

ggplot(combined_data, aes(ageofonset/12, fill=Dataset)) + 
  geom_histogram(binwidth = 2, color="black", position="dodge", alpha=0.7) + 
  geom_vline(xintercept=72/12, linetype='dashed', col = 'red', size=1) +
  theme_minimal() + 
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    legend.title = element_text(size=12),
    legend.text = element_text(size=11),
    plot.title = element_text(size=16, face="bold", hjust=0.5)
  ) +
  scale_x_continuous(breaks = seq(0, max(combined_data$ageofonset/12, na.rm=TRUE), by=10)) +
  scale_fill_manual(values = c(
    "Training" = "#93c47d",
    "Validation" = "#6d9eeb",
    "Test" = "#ea9999",
    "WT Control" = "#f6b26b"
  )) +
  labs(
    title = "Age Distribution of Cancer Onset",
    y = "Number of Individuals", 
    x = "Age of First Cancer Onset (years)",
    fill = "Cohort"
  )


#######################################################
######### Figure 2C (Age Sample vs AgeOnset) ##########
#######################################################

ggplot(combined_data, aes(ageofonset/12, agesamplecollection/12, color=Dataset)) + 
  geom_point() + 
  geom_vline(xintercept=72/12, linetype='dashed', col = 'red', size=1) +
  theme_minimal() + 
  scale_color_manual(values = c(
    "Training" = "#93c47d",
    "Validation" = "#6d9eeb",
    "Test" = "#ea9999",
    "WT Control" = "#f6b26b"
  )) +
  labs(
    y = "Age of Sample Collection (years)", 
    x = "Age of First Cancer Onset (years)",
    color = "Cohort"
  )
