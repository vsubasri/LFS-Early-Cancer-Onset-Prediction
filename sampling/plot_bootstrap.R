library(ggplot2)
library(dplyr)

dir <- "/Users/vallijahsubasri/Documents/lfs_ageofonset/predictions/"

## Helper function to calculate NPV
calculate_npv <- function(predictions, labels, cutoff) {
  pred_neg <- predictions <= cutoff  # predicted negatives
  actual_neg <- labels == 0  # actual negatives
  tn <- sum(pred_neg & actual_neg)  # true negatives
  fn <- sum(!pred_neg & labels == 1)  # false negatives
  return(tn / (tn + fn))
}

## Helper function to process ROC objects
process_roc <- function(results, pred_col, cutoff, model_name, model_type) {
  # Remove duplicates and ensure binary classification
  results <- results[!duplicated(results$ids),]
  
  # Convert labels to binary format and ensure numeric
  results$test_label <- as.numeric(as.factor(results$test_label)) - 1
  
  # Ensure predictions are numeric probabilities
  pred_col_full <- paste0(pred_col, ".Yes")
  if (!pred_col_full %in% names(results)) {
    stop(sprintf("Prediction column %s not found", pred_col_full))
  }
  
  # Add error checking for class counts
  label_counts <- table(results$test_label)
  if (length(label_counts) != 2) {
    warning(sprintf("Found %d classes instead of 2: %s", 
                   length(label_counts), 
                   paste(names(label_counts), collapse=", ")))
    return(list(
      ROCobj = NULL,
      auprc = NULL,
      auc_ci = c(NA, NA),
      f1_ci = c(NA, NA),
      npv_ci = c(NA, NA),
      sens_ci = c(NA, NA),
      spec_ci = c(NA, NA)
    ))
  }
  
  tryCatch({
    # Calculate ROC with confidence intervals
    ROCobj <- if(is.null(cutoff)) {
      ROCInfo_atopt_ci(results, pred_col, "test_label", 1, 1, model_name, model_type)
    } else {
      ROCInfo_atcutoff_ci(results, pred_col, "test_label", cutoff, model_name, model_type)
    }
    
    # Calculate AUPRC
    auprc <- tryCatch({
      pr.curve(
        scores.class0 = results[[pred_col_full]][results$test_label == 0],
        scores.class1 = results[[pred_col_full]][results$test_label == 1]
      )
    }, error = function(e) {
      warning(sprintf("AUPRC calculation failed: %s", e$message))
      return(NULL)
    })
    
    return(list(
      ROCobj = ROCobj,
      auprc = auprc,
      auc_ci = if(!is.null(ROCobj)) ROCobj$auc_ci else c(NA, NA),
      f1_ci = if(!is.null(ROCobj)) ROCobj$f1_ci else c(NA, NA),
      npv_ci = if(!is.null(ROCobj)) ROCobj$npv_ci else c(NA, NA),
      sens_ci = if(!is.null(ROCobj)) ROCobj$sens_ci else c(NA, NA),
      spec_ci = if(!is.null(ROCobj)) ROCobj$spec_ci else c(NA, NA)
    ))
    
  }, error = function(e) {
    warning(sprintf("ROC calculation failed: %s", e$message))
    return(list(
      ROCobj = NULL,
      auprc = NULL,
      auc_ci = c(NA, NA),
      f1_ci = c(NA, NA),
      npv_ci = c(NA, NA),
      sens_ci = c(NA, NA),
      spec_ci = c(NA, NA)
    ))
  })
}

## Function to get all AUC metrics
get_all_auc <- function(dir, n_bootstrap = NULL) {
  # Ensure consistent metric naming
  metrics <- c("auc", "sens", "spec", "f1", "auprc", "npv")
  metric_indices <- c(auc = 6, sens = 7, spec = 8, f1 = 9, auprc = 10, npv = NA)
  test_metrics <- setNames(vector("list", length(metrics)), metrics)
  val_metrics <- setNames(vector("list", length(metrics)), metrics)
  
  # Initialize CI storage with consistent naming for both val and test
  ci_metrics <- list(
    auc_val_ci = list(),
    f1_val_ci = list(),
    npv_val_ci = list(),
    sens_val_ci = list(),
    spec_val_ci = list(),
    auc_test_ci = list(),
    f1_test_ci = list(),
    npv_test_ci = list(),
    sens_test_ci = list(),
    spec_test_ci = list()
  )
  
  files <- list.files(dir, pattern = '_ageofonset_results.rds')
  if (!is.null(n_bootstrap)) {
    files <- head(files, n_bootstrap)
  }
  
  for (file in files) {
    model <- str_replace_all(file, 'NoobCorrected_beta_ProjPC2Adj_|_svmRadial_ageofonset_results.rds', '')
    cat(paste0("Processing model: ", model, '\n'))
    name <- str_replace(file, 'svmRadial_ageofonset_results.rds', 'FNW1_svmRadial_ageofonset_')
    model_name <- unlist(str_split(name, "_"))[length(unlist(str_split(name, "_"))) - 2]
    
    if (file.exists(paste0(dir, name, 'ROCInfoVal.rds'))) {
      # Process validation results with error handling
      val_results <- tryCatch({
        readRDS(paste0(dir, name, 'ROCInfoVal.rds'))[[1]]
      }, error = function(e) NULL)
      
      test_results <- tryCatch({
        readRDS(paste0(dir, name, 'ROCInfoTest.rds'))[[1]]
      }, error = function(e) NULL)
      
      if(is.null(val_results) || is.null(test_results)) {
        warning(sprintf("Could not read results for model: %s", model))
        next
      }
      
      # Get optimal cutoff with error handling
      cutoff <- tryCatch({
        ROCobj_val[[5]]
      }, error = function(e) 0.5)
      if(is.null(cutoff) || is.na(cutoff)) cutoff <- 0.5
      
      # Process both sets
      val_processed <- process_roc(val_results, 'test_pred_calibrated', cutoff, model_name, model_type)
      test_processed <- process_roc(test_results, 'test_pred_calibrated', cutoff, model_name, model_type)
      
      # Extract metrics with safe access
      for (metric in metrics) {
        if (!is.na(metric_indices[metric])) {
          idx <- metric_indices[metric]
          # Handle AUC/AUROC specially with error checking
          if (metric == "auc") {
            test_metrics[[metric]][model] <- if(!is.null(test_processed$ROCobj) && length(test_processed$ROCobj) >= idx) {
              test_processed$ROCobj[[idx]]
            } else NA
            val_metrics[[metric]][model] <- if(!is.null(val_processed$ROCobj) && length(val_processed$ROCobj) >= idx) {
              val_processed$ROCobj[[idx]]
            } else NA
          } else {
            # Handle other metrics
            test_metrics[[metric]][model] <- if(!is.null(test_processed$ROCobj) && length(test_processed$ROCobj) >= idx) {
              test_processed$ROCobj[[idx]]
            } else NA
            val_metrics[[metric]][model] <- if(!is.null(val_processed$ROCobj) && length(val_processed$ROCobj) >= idx) {
              val_processed$ROCobj[[idx]]
            } else NA
          }
        }
      }
      
      # Handle AUPRC separately
      test_metrics$auprc[model] <- if(!is.null(test_processed$auprc)) {
        test_processed$auprc$auc.integral
      } else NA
      val_metrics$auprc[model] <- if(!is.null(val_processed$auprc)) {
        val_processed$auprc$auc.integral
      } else NA
      
      # Calculate NPV with error handling
      test_metrics$npv[model] <- tryCatch({
        calculate_npv(test_results$test_pred_calibrated.Yes, test_results$test_label, cutoff)
      }, error = function(e) NA)
      val_metrics$npv[model] <- tryCatch({
        calculate_npv(val_results$test_pred_calibrated.Yes, val_results$test_label, cutoff)
      }, error = function(e) NA)
      
      # Store all confidence intervals for both validation and test
      ci_metrics$auc_val_ci[[model]] <- val_processed$auc_ci
      ci_metrics$f1_val_ci[[model]] <- val_processed$f1_ci
      ci_metrics$npv_val_ci[[model]] <- val_processed$npv_ci
      ci_metrics$sens_val_ci[[model]] <- val_processed$sens_ci
      ci_metrics$spec_val_ci[[model]] <- val_processed$spec_ci
      
      ci_metrics$auc_test_ci[[model]] <- test_processed$auc_ci
      ci_metrics$f1_test_ci[[model]] <- test_processed$f1_ci
      ci_metrics$npv_test_ci[[model]] <- test_processed$npv_ci
      ci_metrics$sens_test_ci[[model]] <- test_processed$sens_ci
      ci_metrics$spec_test_ci[[model]] <- test_processed$spec_ci
    }
  }
  
  # Initialize results with correct column names including both val and test CIs
  results <- data.frame(matrix(NA, nrow=length(files), ncol=24))
  colnames(results) <- c(
    'auc_val', 'sens_val', 'spec_val', 'f1_val', 'auprc_val', 'npv_val',
    'auc_test', 'sens_test', 'spec_test', 'f1_test', 'auprc_test', 'npv_test',
    'auc_val_ci_lower', 'auc_val_ci_upper',
    'f1_val_ci_lower', 'f1_val_ci_upper',
    'npv_val_ci_lower', 'npv_val_ci_upper',
    'auc_test_ci_lower', 'auc_test_ci_upper',
    'f1_test_ci_lower', 'f1_test_ci_upper',
    'npv_test_ci_lower', 'npv_test_ci_upper'
  )
  
  # Process metrics
  for (metric in metrics) {
    val_col <- paste0(metric, "_val")
    test_col <- paste0(metric, "_test")
    results[[val_col]] <- unlist(val_metrics[[metric]])
    results[[test_col]] <- unlist(test_metrics[[metric]])
  }
  
  # Process CIs for both validation and test
  results$auc_val_ci_lower <- sapply(ci_metrics$auc_val_ci, function(x) if(length(x) >= 1) x[1] else NA)
  results$auc_val_ci_upper <- sapply(ci_metrics$auc_val_ci, function(x) if(length(x) >= 2) x[2] else NA)
  results$f1_val_ci_lower <- sapply(ci_metrics$f1_val_ci, function(x) if(length(x) >= 1) x[1] else NA)
  results$f1_val_ci_upper <- sapply(ci_metrics$f1_val_ci, function(x) if(length(x) >= 2) x[2] else NA)
  results$npv_val_ci_lower <- sapply(ci_metrics$npv_val_ci, function(x) if(length(x) >= 1) x[1] else NA)
  results$npv_val_ci_upper <- sapply(ci_metrics$npv_val_ci, function(x) if(length(x) >= 2) x[2] else NA)
  
  results$auc_test_ci_lower <- sapply(ci_metrics$auc_test_ci, function(x) if(length(x) >= 1) x[1] else NA)
  results$auc_test_ci_upper <- sapply(ci_metrics$auc_test_ci, function(x) if(length(x) >= 2) x[2] else NA)
  results$f1_test_ci_lower <- sapply(ci_metrics$f1_test_ci, function(x) if(length(x) >= 1) x[1] else NA)
  results$f1_test_ci_upper <- sapply(ci_metrics$f1_test_ci, function(x) if(length(x) >= 2) x[2] else NA)
  results$npv_test_ci_lower <- sapply(ci_metrics$npv_test_ci, function(x) if(length(x) >= 1) x[1] else NA)
  results$npv_test_ci_upper <- sapply(ci_metrics$npv_test_ci, function(x) if(length(x) >= 2) x[2] else NA)
  
  return(results)
}

aggregate_bootstrap <- function(bsdir, location, final_ROCobj_val, final_ROCobj_test, n_bootstrap = NULL) {
    # For random sampling models, use existing logic
    bootstrap <- get_all_auc(bsdir, n_bootstrap)
    
    # Get metrics directly from ROC objects
    final_metrics <- c(
        auc_val = final_ROCobj_val[[6]],
        sens_val = final_ROCobj_val[[7]],
        spec_val = final_ROCobj_val[[8]],
        f1_val = final_ROCobj_val[[9]],
        auprc_val = final_ROCobj_val[[10]],
        npv_val = final_ROCobj_val$npv,
        auc_test = final_ROCobj_test[[6]],
        sens_test = final_ROCobj_test[[7]],
        spec_test = final_ROCobj_test[[8]],
        f1_test = final_ROCobj_test[[9]],
        auprc_test = final_ROCobj_test[[10]],
        npv_test = final_ROCobj_test$npv,
        # Validation set CIs
        auc_val_ci_lower = final_ROCobj_val$auc_ci[1],
        auc_val_ci_upper = final_ROCobj_val$auc_ci[2],
        f1_val_ci_lower = final_ROCobj_val$f1_ci[1],
        f1_val_ci_upper = final_ROCobj_val$f1_ci[2],
        npv_val_ci_lower = final_ROCobj_val$npv_ci[1],
        npv_val_ci_upper = final_ROCobj_val$npv_ci[2],
        # Test set CIs
        auc_test_ci_lower = final_ROCobj_test$auc_ci[1],
        auc_test_ci_upper = final_ROCobj_test$auc_ci[2],
        f1_test_ci_lower = final_ROCobj_test$f1_ci[1],
        f1_test_ci_upper = final_ROCobj_test$f1_ci[2],
        npv_test_ci_lower = final_ROCobj_test$npv_ci[1],
        npv_test_ci_upper = final_ROCobj_test$npv_ci[2]
    )
    
    bootstrap <- rbind(bootstrap, final_metrics)
    bootstrap$type <- c(rep("random", nrow(bootstrap) - 1), "lfs")
    bootstrap$func <- location
    
    return(bootstrap)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Modify the bootstrap combination code to handle NULL returns
bootstrap_dir <- "/Users/vallijahsubasri/Documents/lfs_ageofonset/random_sampling/"
n_bootstrap <- 100  # or whatever number you want to use

bs_results <- list(
  tss200 = aggregate_bootstrap(paste0(bootstrap_dir, "bootstrap_tss200/"), "TSS200", final_ROCobj_val, final_ROCobj_test, n_bootstrap),
  body = aggregate_bootstrap(paste0(bootstrap_dir, "bootstrap_body/"), "Body", final_ROCobj_val, final_ROCobj_test, n_bootstrap),
  utr3 = aggregate_bootstrap(paste0(bootstrap_dir, "bootstrap_3utr/"), "3UTR", final_ROCobj_val, final_ROCobj_test, n_bootstrap),
  utr5 = aggregate_bootstrap(paste0(bootstrap_dir, "bootstrap_5utr/"), "5UTR", final_ROCobj_val, final_ROCobj_test, n_bootstrap),
  tss1500 = aggregate_bootstrap(paste0(bootstrap_dir, "bootstrap_tss1500/"), "TSS1500", final_ROCobj_val, final_ROCobj_test, n_bootstrap),
  exon1 = aggregate_bootstrap(paste0(bootstrap_dir, "bootstrap_1stexon/"), "1stExon", final_ROCobj_val, final_ROCobj_test, n_bootstrap),
  lfs = aggregate_bootstrap(paste0(bootstrap_dir, "bootstrap_lfs/"), "lfs", final_ROCobj_val, final_ROCobj_test, n_bootstrap)
)

# Combine results, filtering out NULLs
bootstrap <- do.call(rbind, Filter(Negate(is.null), bs_results))

# Define metrics and process data
metric_names <- c(
  "auc_val", "sens_val", "spec_val", "f1_val", "auprc_val", "npv_val",
  "auc_test", "sens_test", "spec_test", "f1_test", "auprc_test", "npv_test",
  "auc_val_ci_lower", "auc_val_ci_upper", "f1_val_ci_lower", "f1_val_ci_upper", "npv_val_ci_lower", "npv_val_ci_upper",
  "auc_test_ci_lower", "auc_test_ci_upper", "f1_test_ci_lower", "f1_test_ci_upper", "npv_test_ci_lower", "npv_test_ci_upper",
  "type", "func"
)
colnames(bootstrap) <- metric_names

# Process the data without the problematic NPV handling
bootstrap_processed <- bootstrap %>%
  mutate(across(matches("ci_"), as.numeric)) %>%
  reshape2::melt(id.vars = c("type", "func", 
                          # Include both validation and test CIs
                          "auc_val_ci_lower", "auc_val_ci_upper",
                          "f1_val_ci_lower", "f1_val_ci_upper",
                          "npv_val_ci_lower", "npv_val_ci_upper",
                          "auc_test_ci_lower", "auc_test_ci_upper",
                          "f1_test_ci_lower", "f1_test_ci_upper",
                          "npv_test_ci_lower", "npv_test_ci_upper"),
                variable.name = "metric",
                value.name = "value") %>%
  mutate(
    dataset = factor(ifelse(grepl("_val$", metric), "Validation", "Test"),
                    levels = c("Validation", "Test")),
    metric_clean = factor(case_when(
      grepl("^auc", metric) ~ "AUROC",
      grepl("^f1", metric) ~ "F1",
      grepl("^npv", metric) ~ "NPV",
      TRUE ~ NA_character_
    ), levels = c("AUROC", "F1", "NPV")),
    ci_lower = case_when(
      metric_clean == "AUROC" & dataset == "Validation" ~ auc_val_ci_lower,
      metric_clean == "AUROC" & dataset == "Test" ~ auc_test_ci_lower,
      metric_clean == "F1" & dataset == "Validation" ~ f1_val_ci_lower,
      metric_clean == "F1" & dataset == "Test" ~ f1_test_ci_lower,
      metric_clean == "NPV" & dataset == "Validation" ~ npv_val_ci_lower,
      metric_clean == "NPV" & dataset == "Test" ~ npv_test_ci_lower
    ),
    ci_upper = case_when(
      metric_clean == "AUROC" & dataset == "Validation" ~ auc_val_ci_upper,
      metric_clean == "AUROC" & dataset == "Test" ~ auc_test_ci_upper,
      metric_clean == "F1" & dataset == "Validation" ~ f1_val_ci_upper,
      metric_clean == "F1" & dataset == "Test" ~ f1_test_ci_upper,
      metric_clean == "NPV" & dataset == "Validation" ~ npv_val_ci_upper,
      metric_clean == "NPV" & dataset == "Test" ~ npv_test_ci_upper
    ),
    type = factor(type, levels = c("lfs", "random"), 
                 labels = c("LFS", "Random")),
    func = factor(func, levels = c("TSS200", "TSS1500", "5UTR", 
                                 "1stExon", "Body", "3UTR", "lfs")),
    value = pmin(pmax(as.numeric(value), 0), 1)  # Clamp values between 0 and 1
  ) %>%
  filter(!is.na(metric_clean)) %>%
  group_by(dataset, metric_clean, type, func) %>%
  mutate(
    ci_lower = pmax(0, ci_lower),
    ci_upper = pmin(1, ci_upper)
  ) %>%
  ungroup()

# Handle NPV values - single step
bootstrap_processed <- bootstrap_processed %>%
  group_by(dataset, func) %>%
  mutate(
    value = if_else(
      metric_clean == "NPV" & type == "LFS",
      max(value[type == "LFS" & metric_clean == "NPV"], na.rm = TRUE),
      value
    )
  ) %>%
  ungroup()

# Create plot with properly aligned elements
dodge_width <- 0.7  # Width for dodging between groups
box_width <- 0.8    # Width of the boxplots

# Add horizontal offsets for each function within groups
bootstrap_processed <- bootstrap_processed %>%
  group_by(metric_clean, dataset) %>%
  mutate(
    x_offset = (as.numeric(func) - mean(as.numeric(func))) * 0.1,
    x_position = as.numeric(metric_clean) + x_offset
  ) %>%
  ungroup()

# Create single position dodge object
pos_dodge <- position_dodge(width = dodge_width)

# Prepare data for error bars with same offset
error_bar_data <- bootstrap_processed %>%
  filter(!is.na(ci_lower) & !is.na(ci_upper))

plot_metrics <- ggplot() +
  facet_wrap(~dataset) +
  # Boxplots for Random data (colored by func)
  geom_boxplot(
    data = subset(bootstrap_processed, type == "Random"),
    aes(x = metric_clean, y = value, fill = func, 
        group = interaction(metric_clean, func)),
    position = pos_dodge,
    width = box_width,
    alpha = 0.6,
    outlier.shape = NA
  ) +
  # Points with manual offset (colored by type)
  geom_point(
    data = bootstrap_processed,
    aes(x = x_position,
        y = value,
        color = type,    # Color by type
        size = type),
    position = position_identity()
  ) +
  # Error bars with matching position (colored by type)
  geom_errorbar(
    data = error_bar_data,
    aes(x = x_position,
        y = value,
        ymin = ci_lower,
        ymax = ci_upper,
        color = type),   # Color by type to match points
    position = position_identity(),
    width = 0.1
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  # Use different color scales for fill (boxplots) and color (points/errors)
  scale_fill_brewer(palette = "Set2") +
  scale_color_manual(values = c("Random" = "darkgray", "LFS" = "black")) +
  scale_shape_manual(values = c("Random" = 16, "LFS" = 18)) +
  scale_size_manual(values = c("Random" = 1.5, "LFS" = 1.5)) +
  guides(
    fill = guide_legend(title = "Genomic Region", order = 1),
    color = guide_legend(title = "Type", order = 2),
    shape = guide_legend(title = "Type", order = 2),
    size = guide_legend(title = "Type", order = 2)
  ) +
  labs(
    y = "Performance",
    x = "Metric"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(plot_metrics)


