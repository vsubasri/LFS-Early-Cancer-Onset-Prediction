library(wesanderson)
library(ggplot2)
library(ggsci)
library(reshape2)

dir <- '/Users/vallijahsubasri/Documents/lfs_ageofonset/predictions/'
setwd(dir)

models <- c("enet","rf","gbm","xgboost","svmLinear","svmRadial","svmPoly","nnet")
models_results <- list()

for (model_type in models) {
  print(model_type)
  ROCobj_val <- readRDS(paste0(dir,'NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_72_FNW1_',model_type,'_ageofonset_ROCInfoVal.rds'))
  val_results <- ROCobj_val[[1]]
  val_results <- val_results[!duplicated(val_results$ids),]
  weight=1
  pred_col='test_pred_calibrated'
  
  # Test set evaluation
  ROCobj_test <- readRDS(paste0(dir,'NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_72_FNW1_',model_type,'_ageofonset_ROCInfoTest.rds')) 
  test_results <- ROCobj_test[[1]]
  test_results <- test_results[!duplicated(test_results$ids),]
  
  # Ensure we're using the correct prediction column
  if (model_type %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet","svmPoly")) {
    pred_col <- paste0(pred_col,".Yes")
  }
  
  # Calculate confusion matrix for NPV
  predicted_class <- ifelse(test_results[[pred_col]] >= cutoff, 1, 0)
  actual_class <- test_results$test_label
  
  # Create confusion matrix
  cm <- table(actual_class, predicted_class)
  print("Confusion Matrix:")
  print(cm)
  
  # Calculate metrics
  tn <- cm[1,1]  # True Negatives
  fn <- cm[2,1]  # False Negatives
  tp <- cm[2,2]  # True Positives
  fp <- cm[1,2]  # False Positives
  
  # Calculate all metrics
  npv <- if((tn + fn) > 0) tn / (tn + fn) else NA
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  precision <- tp / (tp + fp)
  f1 <- if(sensitivity + precision > 0) 2 * (sensitivity * precision) / (sensitivity + precision) else 0
  auc <- ROCobj_test$auc
  
  print(paste("NPV:", npv))
  print(paste("Sensitivity:", sensitivity))
  print(paste("Specificity:", specificity))
  print(paste("F1:", f1))
  print(paste("AUC:", auc))
  
  # Store results
  test = data.frame(
    auc = auc,
    sens = sensitivity,
    spec = specificity,
    f1 = f1,
    npv = npv
  )
  
  # Round after all calculations are done
  test[] <- lapply(test, round, digits=2)
  
  rownames(test) <- c("test")
  models_results[[model_type]] <- test
}

# Combine results
models_results <- do.call(rbind, models_results)
models_results$model <- rownames(models_results)

# Create empty lists to store bootstrap results
bootstrap_results <- list()
n_bootstrap <- 100  # Number of bootstrap iterations

# Bootstrap for standard errors
for (model_type in models) {
  print(model_type)
  ROCobj_test <- readRDS(paste0(dir,'NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_72_FNW1_',model_type,'_ageofonset_ROCInfoTest.rds')) 
  test_results <- ROCobj_test[[1]]
  test_results <- test_results[!duplicated(test_results$ids),]
  
  # Ensure correct prediction column
  if (model_type %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet","svmPoly")) {
    pred_col <- paste0("test_pred_calibrated",".Yes")
  } else {
    pred_col <- "test_pred_calibrated"
  }
  
  boot_metrics <- matrix(NA, nrow=n_bootstrap, ncol=5)
  colnames(boot_metrics) <- c("auc", "sens", "spec", "f1", "npv")
  
  for(i in 1:n_bootstrap) {
    # Sample with replacement
    boot_indices <- sample(1:nrow(test_results), replace=TRUE)
    boot_data <- test_results[boot_indices,]
    
    # Calculate AUC for this bootstrap sample
    pred <- prediction(boot_data[[pred_col]], boot_data$test_label)
    perf <- performance(pred, "auc")
    boot_auc <- as.numeric(perf@y.values[[1]])
    
    # Calculate other metrics
    predicted_class <- ifelse(boot_data[[pred_col]] >= cutoff, 1, 0)
    actual_class <- boot_data$test_label
    cm <- table(actual_class, predicted_class)
    
    if(nrow(cm) < 2 || ncol(cm) < 2) {
      next
    }
    
    tn <- cm[1,1]
    fn <- cm[2,1]
    tp <- cm[2,2]
    fp <- cm[1,2]
    
    npv <- if((tn + fn) > 0) tn / (tn + fn) else NA
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    precision <- tp / (tp + fp)
    f1 <- if(sensitivity + precision > 0) 2 * (sensitivity * precision) / (sensitivity + precision) else 0
    
    boot_metrics[i,] <- c(boot_auc, sensitivity, specificity, f1, npv)
  }
  
  # Calculate standard errors
  se <- apply(boot_metrics, 2, sd, na.rm=TRUE)
  bootstrap_results[[model_type]] <- se
}

# Convert bootstrap results to data frame
se_results <- do.call(rbind, bootstrap_results)

# Prepare data for plotting
models_results_plot <- melt(models_results[, c("model", "auc", "sens", "spec", "f1", "npv")],
                          id.vars = "model")

# Convert standard errors to long format
se_long <- melt(se_results)
colnames(se_long) <- c("model", "variable", "se")

# Merge metrics with standard errors
models_results_plot <- merge(models_results_plot, se_long,
                           by = c("model", "variable"))

# Color palette
cbPalette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
               "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")

# Create plot
p <- ggplot(models_results_plot, aes(variable, value, fill=model)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), color="black") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                position=position_dodge(width=0.9),
                width=0.25) +
  scale_fill_manual(values=cbPalette, name="Model") +
  theme_minimal() +
  xlab("Performance Metric") +
  ylab("Value") +
  scale_x_discrete(labels = c("auc" = "AUROC",
                             "sens" = "Sensitivity",
                             "spec" = "Specificity",
                             "f1" = "F1-score",
                             "npv" = "NPV")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display plot
print(p)

