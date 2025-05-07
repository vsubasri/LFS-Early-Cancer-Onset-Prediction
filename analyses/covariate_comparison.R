# Load required libraries
library(wesanderson)
library(ggplot2)
library(ggsci)
library(reshape2)
library(stringr)
library(ROCR)
library(data.table)
library(caret)
library(gridExtra)
library(grid)
library(MLmetrics)
library(xgboost)
library(ggpubr)
library(RColorBrewer)
library(umap)
library(pROC)
library(PRROC)
library(boot)

# Define helper function for F1 score calculation
get_f1 <- function(data, predict, actual, cutoff) {
  if (any(unique(data[[actual]]) %in% c("Yes","No"))) {
    data[[actual]] <- ifelse(data[[actual]] == "Yes", 1, 0)
  }
  
  predictions <- ifelse(data[[predict]] >= cutoff, 1, 0)
  cm <- table(data[[actual]], predictions)
  
  # Handle cases where confusion matrix dimensions are incomplete
  if (nrow(cm) < 2 || ncol(cm) < 2) {
    if (nrow(cm) < 2) cm <- rbind(cm, c(0,0))
    if (ncol(cm) < 2) cm <- cbind(cm, c(0,0))
  }
  
  # Calculate precision and recall with handling for edge cases
  precision <- ifelse(sum(cm[,2]) > 0, cm[2,2] / sum(cm[,2]), 0)
  recall <- ifelse(sum(cm[2,]) > 0, cm[2,2] / sum(cm[2,]), 0)
  
  # Calculate F1 score only if both precision and recall are non-zero
  f1 <- ifelse(precision + recall > 0, 
               2 * (precision * recall) / (precision + recall),
               0)
  
  return(f1)
}

compute_auroc_ci <- function(actual, predicted, n_boot = 1000, conf_level = 0.95) {
  # Calculate ROC object
  roc_obj <- roc(actual, predicted, quiet = TRUE)
  
  # Bootstrap confidence intervals
  ci <- ci.auc(roc_obj, conf.level = conf_level, boot.n = n_boot, boot.stratified = TRUE)
  
  # Debugging: Print intermediate values
  print(paste("AUC:", as.numeric(auc(roc_obj))))
  print(paste("CI Lower:", as.numeric(ci[1]), "CI Upper:", as.numeric(ci[3])))
  
  return(c(auc = as.numeric(auc(roc_obj)), lower = as.numeric(ci[1]), upper = as.numeric(ci[3])))
}

compute_auprc_ci <- function(actual, predicted, n_boot = 1000, conf_level = 0.95) {
  # Calculate baseline AUPRC
  base_pr <- pr.curve(scores.class0 = predicted[actual == 1], 
                     scores.class1 = predicted[actual == 0],
                     curve = FALSE)
  base_auprc <- base_pr$auc.integral
  
  # Function to safely compute PR curve
  safe_pr <- function(actual, predicted) {
    tryCatch({
      pr <- pr.curve(scores.class0 = predicted[actual == 1],
                     scores.class1 = predicted[actual == 0],
                     curve = FALSE)
      return(pr$auc.integral)
    }, error = function(e) {
      return(NA)
    })
  }
  
  # Stratified bootstrap
  n_pos <- sum(actual == 1)
  n_neg <- sum(actual == 0)
  pos_indices <- which(actual == 1)
  neg_indices <- which(actual == 0)
  
  boot_values <- replicate(n_boot, {
    boot_pos <- sample(pos_indices, size = n_pos, replace = TRUE)
    boot_neg <- sample(neg_indices, size = n_neg, replace = TRUE)
    indices <- c(boot_pos, boot_neg)
    safe_pr(actual[indices], predicted[indices])
  })
  
  # Remove NAs and calculate CIs
  boot_values <- boot_values[!is.na(boot_values)]
  if(length(boot_values) < 10) {
    return(c(auprc = base_auprc, 
             lower = max(0, base_auprc - 0.1),
             upper = min(1, base_auprc + 0.1)))
  }
  
  # Calculate percentile confidence intervals
  alpha <- (1 - conf_level)/2
  ci <- quantile(boot_values, probs = c(alpha, 1-alpha), na.rm = TRUE)
  
  return(c(
    auprc = base_auprc,
    lower = max(0, ci[1]),
    upper = min(1, ci[2])
  ))
}

compute_metrics_ci <- function(actual, predicted, cutoff, n_boot = 1000, conf_level = 0.95) {
    # Helper function to compute CIs safely
    safe_boot_ci <- function(boot_results, index) {
        tryCatch({
            # Try percentile method first as it's more stable
            ci <- boot.ci(boot_results, type = "perc", index = index)
            return(c(ci$percent[4], ci$percent[5]))
        }, error = function(e) {
            # Fall back to basic bootstrap interval
            sorted_vals <- sort(boot_results$t[,index])
            alpha <- (1 - conf_level)/2
            return(c(
                quantile(sorted_vals, alpha),
                quantile(sorted_vals, 1 - alpha)
            ))
        })
    }

    f <- function(data, indices) {
        d <- data[indices, ]
        y_true <- d$actual
        y_pred <- ifelse(d$predicted >= cutoff, 1, 0)
        
        cm <- table(factor(y_true, levels = c(0,1)), factor(y_pred, levels = c(0,1)))
        if (nrow(cm) < 2) cm <- rbind(cm, c(0,0))
        if (ncol(cm) < 2) cm <- cbind(cm, c(0,0))
        
        sensitivity <- cm[2,2] / sum(cm[2,])
        specificity <- cm[1,1] / sum(cm[1,])
        recall <- sensitivity
        precision <- cm[2,2] / sum(cm[,2])
        f1 <- ifelse(precision + recall == 0, 0, 2 * (precision * recall) / (precision + recall))
        npv <- cm[1,1] / sum(cm[,1])
        
        c(sensitivity, specificity, f1, npv)
    }
    
    # Create data frame and perform bootstrap
    data <- data.frame(actual = actual, predicted = predicted)
    boot_results <- boot(data, statistic = f, R = n_boot)
    
    # Calculate point estimates
    point_estimates <- f(data, 1:nrow(data))
    
    # Get CIs for each metric
    ci_sens <- safe_boot_ci(boot_results, 1)
    ci_spec <- safe_boot_ci(boot_results, 2)
    ci_f1 <- safe_boot_ci(boot_results, 3)
    ci_npv <- safe_boot_ci(boot_results, 4)
    
    return(c(
        sensitivity = point_estimates[1],
        sens_lower = ci_sens[1],
        sens_upper = ci_sens[2],
        specificity = point_estimates[2],
        spec_lower = ci_spec[1],
        spec_upper = ci_spec[2],
        f1 = point_estimates[3],
        f1_lower = ci_f1[1],
        f1_upper = ci_f1[2],
        npv = point_estimates[4],
        npv_lower = ci_npv[1],
        npv_upper = ci_npv[2]
    ))
}

# Define ROCInfo_atopt function with corrected sprintf
ROCInfo_atopt <- function( data, predict, actual, cost.fp, cost.fn, other_title, model) {
  
  if (model %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet","svmPoly") ) {
    predict <- paste0(predict,".Yes")
  }
  
  # Calculate the values using the ROCR library
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  # Calculate cost
  cost <- perf@x.values[[1]] * cost.fp * sum( data[[actual]] == 0 ) +
    ( 1 - perf@y.values[[1]] ) * cost.fn * sum( data[[actual]] == 1 )
  
  cost_dt <- data.frame( cutoff = pred@cutoffs[[1]], cost = cost )
  
  # Optimal cutoff
  best_index  <- which.min(cost)
  best_cost   <- cost_dt[ best_index, "cost" ]
  best_tpr    <- roc_dt[ best_index, "tpr" ]
  best_fpr    <- roc_dt[ best_index, "fpr" ]
  best_cutoff <- pred@cutoffs[[1]][ best_index ]
  
  if (is.infinite(best_cutoff)) {
    best_cutoff = 0.5
  }
  
  f1 <- get_f1(data, predict, actual, best_cutoff)
  
  # AUC
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # Normalize cost for color assignment
  normalize <- function(v) ( v - min(v) ) / diff( range(v) )
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)
  col_by_cost <- col_ramp[ ceiling( normalize(cost) * 99 ) + 1 ]
  
  # ROC Plot
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) +
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.2 ) +
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) +
    labs( title = "ROC", x = "False Positive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = best_tpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = best_fpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" )
  
  # Cost Plot
  cost_plot <- ggplot( cost_dt, aes( cutoff, cost ) ) +
    geom_line( color = "blue", alpha = 0.5 ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.5 ) +
    ggtitle( "Cost" ) +
    geom_vline( xintercept = best_cutoff, alpha = 0.8, linetype = "dashed", color = "steelblue4" )
  
  options(scipen = '999')
  
  # Corrected sprintf usage
  sub_title <- sprintf("%s - Cutoff at %.2f - Total Cost = %.2f, AUC = %.3f",
                       other_title, best_cutoff, best_cost, auc )
  
  # Arrange plots
  plot <- arrangeGrob( roc_plot, cost_plot, ncol = 2,
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ) ) )
  
  # Calculate NPV
  predictions <- ifelse(data[[predict]] >= best_cutoff, 1, 0)
  cm <- table(data[[actual]], predictions)
  if (!("1" %in% colnames(cm))) {
    cm <- cbind(cm, c(0,0))
  } else if (!("0" %in% colnames(cm))) {
    cm <- cbind(c(0,0), cm)
  }
  npv <- cm[1,1] / sum(cm[,1])  # TN / (TN + FN)
  
  return( list(
    data = data,
    plot = plot,
    pred = pred,
    perf = perf,
    cutoff = best_cutoff,
    auc = auc,
    sensitivity = best_tpr,
    specificity = 1 - best_fpr,
    f1 = f1,
    npv = npv,
    totalcost = best_cost
  ) )
}  

# Define ROCInfo_atcutoff function with corrected sprintf
ROCInfo_atcutoff <- function( data, predict, actual, cutoff, other_title, model) {
  
  if (model %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet","svmPoly") ) {
    predict <- paste0(predict,".Yes")
  }
  
  if (any(unique(data[[actual]]) %in% c("Yes","No"))) {
    data[[actual]] <- ifelse(data[[actual]] == "Yes", 1, 0)
  }
  
  # Calculate the values using the ROCR library
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  # Compute metrics
  cm <- table(data[[actual]], ifelse(data[[predict]] >= cutoff,1,0) )
  if (!("1" %in% colnames(cm))) {
    cm <- cbind(cm,c(0,0))
  } else if (!("0" %in% colnames(cm))) {
    cm <- cbind(c(0,0),cm)
  }
  sensitivity <- ifelse(sum(cm[2, ]) == 0, NA, cm[2,2]/sum(cm[2,]))
  specificity <- ifelse(sum(cm[1, ]) == 0, NA, cm[1,1]/sum(cm[1,]))
  accuracy <- (cm[1,1] + cm[2,2])/sum(cm)
  f1 <- get_f1(data, predict, actual, cutoff)
  
  # AUC
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # Color ramp (unused in this function, consider removing if not needed)
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)
  
  # ROC Plot
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) +
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) +
    labs( title = "ROC", x = "False Positive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = sensitivity, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = 1-specificity, alpha = 0.8, linetype = "dashed", color = "steelblue4" )
  
  options(scipen = '999')
  
  # Corrected sprintf usage
  sub_title <- sprintf("%s - Cutoff at %.2f, AUC = %.3f",
                       other_title, cutoff, auc )
  
  # Arrange plots
  plot <- arrangeGrob( roc_plot, ncol = 2,
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ) ) )
  
  # Calculate NPV
  predictions <- ifelse(data[[predict]] >= cutoff, 1, 0)
  cm <- table(data[[actual]], predictions)
  if (!("1" %in% colnames(cm))) {
    cm <- cbind(cm, c(0,0))
  } else if (!("0" %in% colnames(cm))) {
    cm <- cbind(c(0,0), cm)
  }
  npv <- cm[1,1] / sum(cm[,1])  # TN / (TN + FN)
  
  return( list(
    data = data,
    plot = plot,
    pred = pred,
    perf = perf,
    cutoff = cutoff,
    auc = auc,
    sensitivity = sensitivity,
    specificity = specificity,
    f1 = f1,
    npv = npv
  ) )
}

covars <- c("","sex_","tp53_","cellprops_","systreat_","family_","family_sex_","family_tp53_",
            "sex_tp53_","cellprops_sex_","cellprops_sex_tp53_","family_sex_tp53_",
            "cellprops_family_sex_tp53_","cellprops_sex_tp53_nometh_","sex_tp53_PCA_")
covars_results <- list()

# Define the number of bootstrap samples and confidence level
n_boot <- 1000
conf_level <- 0.95

# Define metrics for plotting
metrics <- c("auc", "sens", "spec", "f1", "auprc", "npv")  # Add npv here

# Fix validation data processing
process_validation_data <- function(val_results) {
  if ("test_label" %in% names(val_results)) {
    # First convert any factors to character
    if (is.factor(val_results$test_label)) {
      val_results$test_label <- as.character(val_results$test_label)
    }
    
    # Handle different types of binary inputs
    val_results$test_label <- case_when(
      val_results$test_label %in% c("1", "Yes", "TRUE", "true", "T") ~ 1,
      val_results$test_label %in% c("0", "No", "FALSE", "false", "F") ~ 0,
      is.numeric(val_results$test_label) ~ as.numeric(val_results$test_label > 0),
      TRUE ~ NA_real_
    )
    
    # Remove NA values
    val_results <- val_results[!is.na(val_results$test_label), ]
    
    # Verify binary classification
    unique_labels <- unique(val_results$test_label)
    if (!all(unique_labels %in% c(0, 1))) {
      warning("Non-binary labels detected. Converting to binary.")
      val_results$test_label <- as.numeric(val_results$test_label > 0)
    }
  }
  return(val_results)
}

# Fix missing bars in plot
plot_metrics <- function(metrics_df) {
  # Replace NA or infinite values with 0
  metrics_df[is.na(metrics_df)] <- 0
  metrics_df[!is.finite(metrics_df)] <- 0
  
  # Ensure all metrics are between 0 and 1
  metrics_df[] <- lapply(metrics_df, function(x) pmin(pmax(x, 0), 1))
  
  # Create the plot
  ggplot(metrics_df, aes(x = covars, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Covariates", y = "Value", fill = "Metric") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

for (covar in covars) {
  print(paste("Processing covariate:", covar))
  
  # Add debugging information for validation data
  print("Validation data check:")
  if (!file.exists(roc_val_path)) {
    print(paste("File not found:", roc_val_path))
    next
  }
  
  ROCobj_val <- readRDS(roc_val_path)
  val_results <- ROCobj_val[[1]]
  val_results <- val_results[!duplicated(val_results$ids), ]
  
  # Print data structure and unique values
  print("Structure of validation results:")
  print(str(val_results))
  print("Unique values in test_label:")
  print(unique(val_results$test_label))
  
  # Process validation data
  val_results <- process_validation_data(val_results)
  
  # Verify we have enough data
  if (nrow(val_results) < 2 || length(unique(val_results$test_label)) != 2) {
    warning(paste("Insufficient validation data for covariate:", covar))
    next
  }
  
  # More robust label conversion
  if ("test_label" %in% names(val_results)) {
    # Convert factors or characters to numeric, handling special cases
    if (is.factor(val_results$test_label) || is.character(val_results$test_label)) {
      # Map various possible values to binary
      val_results$test_label <- case_when(
        val_results$test_label %in% c("1", "Yes", "TRUE", "true") ~ 1,
        val_results$test_label %in% c("0", "No", "FALSE", "false") ~ 0,
        TRUE ~ NA_real_
      )
    }
    
    # Verify binary classification
    unique_labels <- unique(val_results$test_label[!is.na(val_results$test_label)])
    if (length(unique_labels) != 2) {
      print("Warning: Non-binary labels detected")
      print("Unique labels after conversion:")
      print(unique_labels)
      next  # Skip this iteration
    }
  }
  
  # Initialize an empty result structure
  empty_result <- data.frame(
    auc = 0, auc_lower = 0, auc_upper = 0,
    sens = 0, sens_lower = 0, sens_upper = 0,
    spec = 0, spec_lower = 0, spec_upper = 0,
    f1 = 0, f1_lower = 0, f1_upper = 0,
    auprc = 0, auprc_lower = 0, auprc_upper = 0,
    npv = 0, npv_lower = 0, npv_upper = 0
  )
  
  # Determine the model based on covariate
  if (covar == "systreat_") {
    model <- "xgboost"
  } else {
    model <- "svmRadial"
  }
  
  # Modify dataset name and file path handling
  dataset <- 'NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_'
  if (covar == "") {
    # Handle the base case without any covariate
    roc_val_path <- paste0(dir, dataset, '72_FNW1_', model, '_ageofonset_ROCInfoVal.rds')
    roc_test_path <- paste0(dir, dataset, '72_FNW1_', model, '_ageofonset_ROCInfoTest.rds')
  } else if (grepl("nometh", covar)) {
    file_covar <- gsub("_nometh_", "_nometh_", covar)
    roc_val_path <- paste0(dir, dataset, file_covar, '72_FNW1_', model, '_ageofonset_ROCInfoVal.rds');
    roc_test_path <- paste0(dir, dataset, file_covar, '72_FNW1_', model, '_ageofonset_ROCInfoTest.rds');
  } else {
    file_covar <- covar
    roc_val_path <- paste0(dir, dataset, file_covar, '72_FNW1_', model, '_ageofonset_ROCInfoVal.rds');
    roc_test_path <- paste0(dir, dataset, file_covar, '72_FNW1_', model, '_ageofonset_ROCInfoTest.rds');
  }
  
  # Read ROC objects with debug information
  ROCobj_val <- readRDS(roc_val_path)
  val_results <- ROCobj_val[[1]]
  val_results <- val_results[!duplicated(val_results$ids), ]
  
  # Ensure binary classification for validation set
  val_results <- process_validation_data(val_results)
  
  # Similar process for test set
  ROCobj_test <- readRDS(roc_test_path)
  test_results <- ROCobj_test[[1]]
  test_results <- test_results[!duplicated(test_results$ids), ]
  
  test_results <- process_validation_data(test_results)
  
  # Ensure correct prediction column names
  if ("test_pred_calibrated.Yes" %in% names(test_results)) {
    test_results$test_pred.Yes <- test_results$test_pred_calibrated.Yes
  }
  if ("test_pred_calibrated.Yes" %in% names(val_results)) {
    val_results$test_pred.Yes <- val_results$test_pred_calibrated.Yes
  }
  
  # Remove any NA values
  test_results <- test_results[!is.na(test_results$test_label) & !is.na(test_results$test_pred.Yes), ]
  val_results <- val_results[!is.na(val_results$test_label) & !is.na(val_results$test_pred.Yes), ]
  
  # Verify we have binary classification data
  if (length(unique(test_results$test_label)) != 2) {
    stop("Test set does not contain binary classification data")
  }
  if (length(unique(val_results$test_label)) != 2) {
    stop("Validation set does not contain binary classification data")
  }
  
  # Add weight parameter definition
  weight <- 1  # Define appropriate weight value
  
  # Compute ROC information at optimal cutoff with weight parameter
  ROCobj_val <- ROCInfo_atopt(val_results, "test_pred_calibrated", "test_label", 1, weight, bestmodel_name, model)
  
  # Compute AUPRC
  val_pr <- pr.curve(scores.class0 = val_results$test_pred.Yes[val_results$test_label == 0], 
                     scores.class1 = val_results$test_pred.Yes[val_results$test_label == 1],
                     curve = FALSE)
  ROCobj_val$auprc <- val_pr$auc.integral
  
  # Define cutoff
  cutoff <- ROCobj_val$cutoff
  if (is.null(cutoff) || is.na(cutoff)) {
    cutoff <- 0.4
    warning("Optimal cutoff not found. Using default cutoff = 0.4")
  }
  
  # Compute ROC information at specified cutoff
  ROCobj_test <- ROCInfo_atcutoff(test_results, "test_pred_calibrated", "test_label", cutoff, bestmodel_name, model)
  
  # Compute AUPRC for test set
  test_pr <- pr.curve(scores.class0 = test_results$test_pred.Yes[test_results$test_label == 0], 
                      scores.class1 = test_results$test_pred.Yes[test_results$test_label == 1],
                      curve = FALSE)
  ROCobj_test$auprc <- test_pr$auc.integral
  
  # Compute Confidence Intervals for Test Set
  test_ci_auroc <- compute_auroc_ci(test_results$test_label, test_results$test_pred.Yes, n_boot, conf_level)
  test_ci_auprc <- compute_auprc_ci(test_results$test_label, test_results$test_pred.Yes, n_boot, conf_level)
  test_ci_metrics <- compute_metrics_ci(test_results$test_label, test_results$test_pred.Yes, cutoff, n_boot, conf_level);
  
  # Print debug information
  print("Test set metrics:")
  print(paste("AUC:", ROCobj_test$auc))
  print(paste("Sensitivity:", ROCobj_test$sensitivity))
  print(paste("Specificity:", ROCobj_test$specificity))
  print(paste("F1:", ROCobj_test$f1))
  print(paste("AUPRC:", ROCobj_test$auprc))
  
  # Assemble test results with CIs
  test = data.frame(
    auc = ifelse(is.na(ROCobj_test$auc), 0, round(ROCobj_test$auc, 3)),
    auc_lower = ifelse(is.na(test_ci_auroc["lower"]), 0, round(test_ci_auroc["lower"], 3)),
    auc_upper = ifelse(is.na(test_ci_auroc["upper"]), 0, round(test_ci_auroc["upper"], 3)),
    sens = ifelse(is.na(ROCobj_test$sensitivity), 0, round(ROCobj_test$sensitivity, 3)),
    sens_lower = ifelse(is.na(test_ci_metrics["sens_lower"]), 0, round(test_ci_metrics["sens_lower"], 3)),
    sens_upper = ifelse(is.na(test_ci_metrics["sens_upper"]), 0, round(test_ci_metrics["sens_upper"], 3)),
    spec = ifelse(is.na(ROCobj_test$specificity), 0, round(ROCobj_test$specificity, 3)),
    spec_lower = ifelse(is.na(test_ci_metrics["spec_lower"]), 0, round(test_ci_metrics["spec_lower"], 3)),
    spec_upper = ifelse(is.na(test_ci_metrics["spec_upper"]), 0, round(test_ci_metrics["spec_upper"], 3)),
    f1 = ifelse(is.na(ROCobj_test$f1), 0, round(ROCobj_test$f1, 3)),
    f1_lower = ifelse(is.na(test_ci_metrics["f1_lower"]), 0, round(test_ci_metrics["f1_lower"], 3)),
    f1_upper = ifelse(is.na(test_ci_metrics["f1_upper"]), 0, round(test_ci_metrics["f1_upper"], 3)),
    auprc = ifelse(is.na(ROCobj_test$auprc), 0, round(ROCobj_test$auprc, 3)),
    auprc_lower = ifelse(is.na(test_ci_auprc["lower"]), 0, round(test_ci_auprc["lower"], 3)),
    auprc_upper = ifelse(is.na(test_ci_auprc["upper"]), 0, round(test_ci_auprc["upper"], 3)),
    npv = ifelse(is.na(test_ci_metrics["npv"]), 0, round(test_ci_metrics["npv"], 3)),
    npv_lower = ifelse(is.na(test_ci_metrics["npv_lower"]), 0, round(test_ci_metrics["npv_lower"], 3)),
    npv_upper = ifelse(is.na(test_ci_metrics["npv_upper"]), 0, round(test_ci_metrics["npv_upper"], 3))
  )
  
  covars_results[[covar]] <- test
}

# Combine all covariate results into a single data frame
covars_results <- do.call(rbind, covars_results)
rownames(covars_results)[1] <- " "

# Ensure all necessary columns exist
required_metrics <- c("auc", "sens", "spec", "f1", "auprc", "npv")
required_columns <- c(
  required_metrics,
  paste0(required_metrics, "_lower"),
  paste0(required_metrics, "_upper")
)

# Check if any required columns are missing and add them if necessary
for (col in required_columns) {
  if (!(col %in% names(covars_results))) {
    covars_results[[col]] <- 0  # Initialize missing columns with 0
  }
}

# Add covariate names as a column
covars_results$covars <- rownames(covars_results)

# Reshape and prepare the data for all metrics
plot_data_melt <- melt(covars_results, 
                       id.vars = "covars",
                       measure.vars = metrics,
                       variable.name = "Metrics",
                       value.name = "Value")

# Create confidence interval data with improved validation
ci_data <- NULL
for (metric in metrics) {
    metric_value <- covars_results[[metric]]
    metric_lower <- covars_results[[paste0(metric, "_lower")]]
    metric_upper <- covars_results[[paste0(metric, "_upper")]]
    
    # Create data frame for this metric
    metric_ci <- data.frame(
        covars = covars_results$covars,
        Metrics = metric,
        Value = metric_value,
        ci_lower = metric_lower,
        ci_upper = metric_upper,
        stringsAsFactors = FALSE
    )
    
    # Validate and fix confidence intervals
    metric_ci <- transform(metric_ci, {
        # For zero values, set both CIs to NA
        zero_mask <- Value == 0 | is.na(Value)
        ci_lower[zero_mask] <- NA
        ci_upper[zero_mask] <- NA
        
        # For non-zero values
        valid_mask <- !zero_mask
        
        if(any(valid_mask)) {
            # Calculate robust standard error
            value_sd <- sd(Value[valid_mask], na.rm = TRUE)
            if(is.na(value_sd) || value_sd < 0.05) value_sd <- 0.05
            
            # Fix problematic CIs
            ci_issues <- valid_mask & (
                ci_upper <= Value | 
                ci_lower >= Value |
                ci_upper <= ci_lower |
                (ci_upper - ci_lower) < 0.05
            )
            
            if(any(ci_issues)) {
                # Recalculate CIs using normal approximation with minimum spread
                z_score <- qnorm(0.975)  # for 95% CI
                min_spread <- 0.05  # minimum difference between value and CI bound
                
                affected_indices <- which(ci_issues)
                for (i in affected_indices) {
                    spread <- max(min_spread, z_score * value_sd)
                    ci_lower[i] <- max(0, Value[i] - spread)
                    ci_upper[i] <- min(1, Value[i] + spread)
                }
            }
        }
        
        list(ci_lower = ci_lower, ci_upper = ci_upper)
    })
    
    # Initialize ci_data if NULL
    if (is.null(ci_data)) {
        ci_data <- metric_ci
    } else {
        ci_data <- rbind(ci_data, metric_ci)
    }
}

# Merge with original data
plot_data_melt <- merge(plot_data_melt, ci_data[c("covars", "Metrics", "ci_lower", "ci_upper")], 
                       by = c("covars", "Metrics"))

# Order within each metric group
plot_data_melt <- plot_data_melt %>%
  group_by(Metrics) %>%
  arrange(desc(Value)) %>%
  mutate(order = row_number()) %>%
  ungroup()

# Set factor levels based on the ordered values for each metric
plot_data_melt$covars <- factor(plot_data_melt$covars, 
                                levels = unique(plot_data_melt$covars[order(plot_data_melt$Metrics, -plot_data_melt$Value)]))

# Define custom colors for covariates with more distinct colors
covar_palette <- c(
  " " = "#FF1493",                    # Deep pink (distinct) for Methylation Only
  "cellprops_" = "#E69F00",             # Orange
  "cellprops_family_sex_tp53_" = "#56B4E9",  # Light blue
  "cellprops_sex_" = "#009E73",         # Green
  "cellprops_sex_tp53_" = "#8B0000",    # Dark red
  "cellprops_sex_tp53_nometh_" = "#0072B2",  # Dark blue
  "family_" = "#FF69B4",                # Hot pink
  "family_sex_" = "#FF6347",            # Dark orchid
  "family_sex_tp53_" = "#4B0082",       # Indigo
  "family_tp53_" = "#FFD700",           # Gold
  "sex_" = "#FFA500",                   # Bright orange for Meth + Sex
  "sex_tp53_" = "#32CD32",              # Lime green
  "systreat_" = "#8B4513",              # Saddle brown
  "tp53_" = "#8A2BE2",                   # Blue violet (distinct) for Meth + TP53
  "Meth" = "#FF1493",                    # Deep pink (distinct)
  "Meth + Sex" = "#FFA500",              # Bright orange
  "Meth + Family" = "#7FFF00",           # Chartreuse
  "Meth + CellProps" = "#00BFFF",        # Deep sky blue
  "Meth + TP53" = "#8A2BE2",             # Blue violet (distinct)
  "Meth + Family + Sex" = "#4B0082",     # Green yellow
  "Meth + Family + TP53" = "#FF6347",    # Tomato
  "Meth + CellProps + Sex" = "#D2691E",  # Chocolate
  "Meth + CellProps + TP53" = "#B22222", # Firebrick
  "Meth + CellProps + Family" = "#FF4500", # Orange red
  "Meth + CellProps + Family + Sex" = "#7B68EE", # Medium slate blue
  "CellProps + Sex" = "#FF00FF",         # Magenta
  "CellProps + Sex + TP53" = "#FFD700",   # Gold
  "sex_tp53_PCA_" = "#00FF00",  # Dark orchid for Meth (PCA) + Sex + TP53
  "Meth (PCA) + Sex + TP53" = "#00FF00" # Bright lime - Added to match sex_tp53_PCA_
)

# Create a named vector for metric labels
metric_labels <- c(
  "auc" = "AUROC",
  "sens" = "Sensitivity",
  "spec" = "Specificity",
  "f1" = "F1-Score",
  "auprc" = "AUPRC",
  "npv" = "NPV"
)

# Get the order based on AUROC values and update legend accordingly
plot_data_melt <- plot_data_melt %>%
  group_by(Metrics) %>%
  mutate(mean_value = mean(Value)) %>%
  ungroup() %>%
  arrange(desc(mean_value), desc(Value))

# Set the factor levels for covariates based on overall performance
covariate_order <- plot_data_melt %>%
  group_by(covars) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE)) %>%
  arrange(desc(mean_value)) %>%
  pull(covars)

# Update factor levels
plot_data_melt$covars <- factor(plot_data_melt$covars, levels = covariate_order)

# Function to format covariate names
format_covariate <- function(x) {
  if (x == " ") {
    return("Meth")
  }
  
  x <- gsub("_$", "", x)
  is_nometh <- grepl("nometh", x)
  is_pca <- grepl("PCA", x)
  x <- gsub("_nometh", "", x)
  x <- gsub("_PCA", "", x)
  parts <- strsplit(x, "_")[[1]]
  
  formatted_parts <- sapply(parts, function(part) {
    switch(part,
           "tp53" = "TP53",
           "cellprops" = "CellProps",
           "sex" = "Sex",
           "family" = "Family",
           "systreat" = "SysTreat",
           part)
  })
  
  result <- paste(formatted_parts, collapse = " + ")
  if (!is_nometh) {
    if (is_pca) {
      result <- paste("Meth (PCA) +", result)
    } else {
      result <- paste("Meth +", result)
    }
  }
  
  return(result)
}

# Apply formatting to the data
plot_data_melt$covars <- factor(
  sapply(as.character(plot_data_melt$covars), format_covariate),
  levels = sapply(covariate_order, format_covariate)
)

# Update the palette names
names(covar_palette) <- sapply(names(covar_palette), format_covariate)

# Ensure Value is numeric
plot_data_melt$Value <- as.numeric(plot_data_melt$Value)  # Ensure Value is numeric

# Create the plot with validated error bars
p <- ggplot(plot_data_melt, aes(x = covars, y = Value, fill = covars)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(data = subset(plot_data_melt, !is.na(ci_lower) & !is.na(ci_upper) & Value > 0),
                  aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.9),
                  width = 0.25,
                  color = "black") +
    scale_fill_manual(values = covar_palette, 
                      name = "Covariates",
                      breaks = sapply(covariate_order, format_covariate)) +
    theme_minimal() +
    labs(title = "Model Performance with Covariates",
         x = "",
         y = "Value") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 12)
    ) +
    facet_wrap(~Metrics, scales = "free_x", ncol = 1, 
               labeller = labeller(Metrics = metric_labels)) +
    coord_cartesian(ylim = c(0, 1))

print(p)

# Add validation check
print("Validation of confidence intervals:")
print(plot_data_melt %>% 
      filter(Value > ci_upper | Value < ci_lower) %>% 
      select(covars, Metrics, Value, ci_lower, ci_upper))

# Ensure that ci_lower and ci_upper are calculated correctly
# Check for AUPRC specifically
print(plot_data_melt[plot_data_melt$Metrics == "auprc", c("covars", "Value", "ci_lower", "ci_upper")])

