library(wesanderson)
library(ggplot2)
library(ggsci)

## get averaged model performance across val and test set
files <- list.files(dir, pattern = "_auc_all.csv", recursive = TRUE, full.names = TRUE)
all <- do.call(rbind, lapply(files, read.csv, row.names = 'X'))
ag <- aggregate(. ~ dataset + preprocessing, all, function(x) c(mean = mean(x,na.rm=TRUE), sd = sd(x,na.rm=TRUE)))
ag_noclin <- ag[!grepl('canceratdraw',ag$preprocessing),]

seed_results <- list()
for (seed in seq(1,5)) {
  print(seed)
  ROCobj_val <- readRDS(paste0(dir,'reps/Seed',seed,'_NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_sex_FNW1_svmRadial_ageofonset_ROCInfoVal.rds'))
  val_results <- ROCobj_val[[1]]
  val_results <- val_results[!duplicated(val_results$ids),]
  weight=1
  pred_col='test_pred_calibrated'
  ROCobj_val <- ROCInfo_atopt(val_results,pred_col,"test_label",1,weight,bestmodel_name, "svmRadial")
  cutoff <- ROCobj_val[[5]]
  val = data.frame(auc = round(ROCobj_val[[6]],2), sens=round(ROCobj_val[[7]],2), spec=round(ROCobj_val[[8]],2),f1=round(ROCobj_val[[9]],2))
  rownames(val) <- c("val")
  ROCobj_test <- readRDS(paste0(dir,'reps/Seed',seed,'_NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_sex_FNW1_svmRadial_ageofonset_ROCInfoTest.rds')) 
  test_results <- ROCobj_test[[1]]
  test_results <- test_results[!duplicated(test_results$ids),]
  ROCobj_test <- ROCInfo_atcutoff(test_results,pred_col,"test_label",cutoff,bestmodel_name,"svmRadial")
  test = data.frame(auc = round(ROCobj_test[[6]],2), sens=round(ROCobj_test[[7]],2), spec=round(ROCobj_test[[8]],2),f1=round(ROCobj_test[[9]],2))
  rownames(test) <- c("test")
  results <- rbind(test)
  seed_results[[seed]] <- results
}
seed_results <- do.call(rbind,seed_results)
seed_results$seed <- rownames(seed_results)
seed_results_plot <- melt(seed_results)

ggplot(seed_results_plot,aes(variable, value, fill=seed)) +
  geom_bar(stat = "identity",position="dodge", colour="black",alpha=0.8) +
  scale_fill_manual(values=wes_palette(n=6, name="IsleofDogs1"),name = "Covariates") + 
  theme_minimal() + 
  xlab("Performance Metric") +
  ylab("Value") 
