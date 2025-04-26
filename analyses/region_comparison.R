library(wesanderson)
library(ggplot2)
library(ggsci)
library(reshape2)

dir <- '/Users/vallijahsubasri/Documents/lfs_ageofonset/predictions/old/'
setwd(dir)

regions <- c("TSS200","TSS1500","Body","gene","3UTR", "5UTR","1stExon")
region_results <- list()

# Function to calculate NPV
calculate_npv <- function(actual, predicted, threshold) {
  # Handle Inf/-Inf values by converting them to numeric
  predicted[is.infinite(predicted)] <- NA
  predicted[is.na(predicted)] <- mean(predicted, na.rm=TRUE)
  
  # Convert predictions to binary using threshold
  pred_binary <- ifelse(predicted >= threshold, 1, 0)
  
  # Create confusion matrix
  tn <- sum(actual == 0 & pred_binary == 0)  # True Negatives
  fn <- sum(actual == 1 & pred_binary == 0)  # False Negatives
  
  # Add detailed debugging
  print("NPV Calculation Details:")
  print(paste("Total samples:", length(actual)))
  print(paste("Threshold:", threshold))
  print(paste("True Negatives:", tn))
  print(paste("False Negatives:", fn))
  print(paste("Total Negative Predictions:", tn + fn))
  
  # Check for zero denominator
  if ((tn + fn) == 0) {
    print("Warning: No negative predictions found")
    return(NA)
  }
  
  # Calculate NPV
  npv <- tn / (tn + fn)
  print(paste("Calculated NPV:", npv))
  return(npv)
}

# Update the ROCInfo_atopt function
ROCInfo_atopt <- function(data, predict, actual, cost.fp, cost.fn, other_title, model) {
  if (model %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet","svmPoly")) {
    predict <- paste0(predict,".Yes")
  }
  
  # Ensure we have valid data
  valid_rows <- !is.na(data[[predict]]) & !is.na(data[[actual]])
  if (sum(valid_rows) == 0) {
    warning("No valid data points found")
    return(list(
      data = data,
      pred = NULL,
      perf = NULL,
      cutoff = 0.5,
      auc = 0,
      sensitivity = 0,
      specificity = 0,
      f1 = 0,
      npv = 0,
      totalcost = 0
    ))
  }
  
  # Use only valid rows for calculations
  clean_data <- data[valid_rows, ]
  
  # Calculate the values using the ROCR library
  pred <- prediction(clean_data[[predict]], clean_data[[actual]])
  perf <- performance(pred, "tpr", "fpr")
  roc_dt <- data.frame(fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
  
  # Calculate cost
  cost <- perf@x.values[[1]] * cost.fp * sum(clean_data[[actual]] == 0) +
    (1 - perf@y.values[[1]]) * cost.fn * sum(clean_data[[actual]] == 1)
  
  cost_dt <- data.frame(cutoff = pred@cutoffs[[1]], cost = cost)
  
  # Optimal cutoff
  best_index <- which.min(cost)
  best_cost <- cost_dt[best_index, "cost"]
  best_tpr <- roc_dt[best_index, "tpr"]
  best_fpr <- roc_dt[best_index, "fpr"]
  best_cutoff <- pred@cutoffs[[1]][best_index]
  
  if (is.infinite(best_cutoff)) {
    best_cutoff = 0.5
  }
  
  # Calculate predictions using clean data
  predictions <- ifelse(clean_data[[predict]] >= best_cutoff, 1, 0)
  
  # Create confusion matrix
  cm <- table(clean_data[[actual]], predictions)
  
  # Handle missing confusion matrix categories
  if (!("1" %in% colnames(cm))) {
    cm <- cbind(cm, c(0,0))
    colnames(cm) <- c("0", "1")
  } else if (!("0" %in% colnames(cm))) {
    cm <- cbind(c(0,0), cm)
    colnames(cm) <- c("0", "1")
  }
  
  # Calculate NPV: TN / (TN + FN)
  npv <- ifelse(sum(cm[,1]) > 0, cm[1,1] / sum(cm[,1]), 0)
  
  # Calculate F1 score
  f1 <- get_f1(clean_data, predict, actual, best_cutoff)
  
  # AUC
  auc <- performance(pred, "auc")@y.values[[1]]
  
  return(list(
    data = data,
    pred = pred,
    perf = perf,
    cutoff = best_cutoff,
    auc = auc,
    sensitivity = best_tpr,
    specificity = 1 - best_fpr,
    f1 = f1,
    npv = npv,
    totalcost = best_cost
  ))
}

# Similarly update ROCInfo_atcutoff function
ROCInfo_atcutoff <- function(data, predict, actual, cutoff, other_title, model) {
  if (model %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet","svmPoly")) {
    predict <- paste0(predict,".Yes")
  }
  
  # Ensure we have valid data
  valid_rows <- !is.na(data[[predict]]) & !is.na(data[[actual]])
  if (sum(valid_rows) == 0) {
    warning("No valid data points found")
    return(list(
      data = data,
      pred = NULL,
      perf = NULL,
      cutoff = cutoff,
      auc = 0,
      sensitivity = 0,
      specificity = 0,
      f1 = 0,
      npv = 0
    ))
  }
  
  # Use only valid rows for calculations
  clean_data <- data[valid_rows, ]
  
  # Calculate predictions using clean data
  predictions <- ifelse(clean_data[[predict]] >= cutoff, 1, 0)
  
  # Create confusion matrix
  cm <- table(clean_data[[actual]], predictions)
  
  # Handle missing confusion matrix categories
  if (!("1" %in% colnames(cm))) {
    cm <- cbind(cm, c(0,0))
    colnames(cm) <- c("0", "1")
  } else if (!("0" %in% colnames(cm))) {
    cm <- cbind(c(0,0), cm)
    colnames(cm) <- c("0", "1")
  }
  
  # Calculate metrics
  sensitivity <- cm[2,2] / sum(cm[2,])
  specificity <- cm[1,1] / sum(cm[1,])
  npv <- ifelse(sum(cm[,1]) > 0, cm[1,1] / sum(cm[,1]), 0)
  
  # Calculate F1 score
  f1 <- get_f1(clean_data, predict, actual, cutoff)
  
  # Calculate AUC
  pred <- prediction(clean_data[[predict]], clean_data[[actual]])
  auc <- performance(pred, "auc")@y.values[[1]]
  
  return(list(
    data = data,
    pred = pred,
    perf = NULL,
    cutoff = cutoff,
    auc = auc,
    sensitivity = sensitivity,
    specificity = specificity,
    f1 = f1,
    npv = npv
  ))
}

for (region in regions) {
  print(region)
  model_type="svmRadial"
  ROCobj_val <- readRDS(paste0(dir,'NoobCorrected_beta_ProjPC2Adj_lfs_',region,'_svmRadial_ageofonset_ROCInfoVal.rds'))
  val_results <- ROCobj_val[[1]]
  val_results <- val_results[!duplicated(val_results$ids),]
  weight=1
  pred_col='test_pred_calibrated'
  ROCobj_val <- ROCInfo_atopt(val_results,pred_col,"test_label",1,weight,bestmodel_name, model_type)
  cutoff <- ROCobj_val$cutoff
  
  ROCobj_test <- readRDS(paste0(dir,'NoobCorrected_beta_ProjPC2Adj_lfs_',region,'_svmRadial_ageofonset_ROCInfoTest.rds')) 
  test_results <- ROCobj_test[[1]]
  test_results <- test_results[!duplicated(test_results$ids),]
  ROCobj_test <- ROCInfo_atcutoff(test_results,pred_col,"test_label",cutoff,bestmodel_name,model_type)
  
  # Create results dataframe with all metrics
  test = data.frame(
    auc = ROCobj_test$auc,
    sens = ROCobj_test$sensitivity,
    spec = ROCobj_test$specificity,
    f1 = ROCobj_test$f1,
    npv = ROCobj_test$npv
  )
  
  # Calculate standard errors
  n <- nrow(test_results)
  test$auc_se <- sqrt(test$auc * (1 - test$auc) / n)
  test$sens_se <- sqrt(test$sens * (1 - test$sens) / n)
  test$spec_se <- sqrt(test$spec * (1 - test$spec) / n)
  test$f1_se <- sqrt(test$f1 * (1 - test$f1) / n)
  test$npv_se <- sqrt(test$npv * (1 - test$npv) / n)
  
  rownames(test) <- c("test")
  region_results[[region]] <- test
}

# Combine results
region_results <- do.call(rbind, region_results)
region_results$region <- rownames(region_results)

# Add debug print to verify values
print("Region Results:")
print(region_results)

# Reshape data for plotting
metrics <- c("auc", "sens", "spec", "f1", "npv")
region_results_plot <- melt(region_results, 
                          measure.vars = metrics,
                          id.vars = "region")

# Fix the standard error addition
se_cols <- paste0(metrics, "_se")
se_data <- melt(region_results[, c("region", se_cols)], 
                id.vars = "region", 
                value.name = "se")
# Match the standard errors with their corresponding metrics
region_results_plot$se <- se_data$se

# Add this debug print before plotting
print("Available metrics in melted data:")
print(unique(region_results_plot$variable))

# Create plot with error bars
ggplot(region_results_plot, aes(x = variable, y = value, fill = region)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                position = position_dodge(width = 0.9),
                width = 0.25) +
  scale_fill_manual(values = pal_jco()(7), name = "Region") + 
  theme_minimal() +
  xlab("Performance Metric") +
  ylab("Value") +
  scale_x_discrete(labels = c("auc" = "AUROC", 
                             "sens" = "Sensitivity",
                             "spec" = "Specificity",
                             "f1" = "F1-score",
                             "npv" = "NPV")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

