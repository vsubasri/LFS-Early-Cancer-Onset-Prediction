library(stringr)
library(ggplot2)
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

get_accuracy <- function(data,predict,actual,cutoff) {
  cm <- data.frame(t(table(data[[actual]], ifelse(data[[predict]] >= cutoff,1,0) )))
  accuracy <- (cm[1,3]+cm[4,3])/sum(cm[,3])
  return(accuracy)
}

opt.cut = function(pred,perf){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

get_auc <- function( data, predict, actual) {
  pred <- prediction( data[[predict]], data[[actual]] )
  auc <- performance( pred, "auc" )@y.values[[1]]
  return(auc)
}

get_f1 <- function (data , predict, actual, cutoff) {
  y_pred <- ifelse(data[[predict]] >= cutoff, 1, 0)
  precision <- Precision(data[[actual]], y_pred, positive = 1)
  recall <- Recall(data[[actual]], y_pred, positive = 1)
  f1 <- 2*precision*recall/(precision+recall)
  return(f1)
}

ROCInfo_atopt_ci <- function( data, predict, actual, cost.fp, cost.fn, other_title, model) {
  
  if (model %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet","svmPoly") ) {
    predict <- paste0(predict,".Yes")
  }
  
  # Validate inputs
  if (!all(data[[actual]] %in% c(0,1,"Yes","No"))) {
    stop("Actual values must be binary (0/1 or Yes/No)")
  }
  
  if (!is.numeric(data[[predict]])) {
    stop("Predictions must be numeric")
  }
  
  # Convert Yes/No to 1/0 if needed
  if (any(data[[actual]] %in% c("Yes","No"))) {
    data[[actual]] <- ifelse(data[[actual]] == "Yes", 1, 0)
  }
  
  # calculate the values using the ROCR library
  # true positive, false postive
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  # cost with the specified false positive and false negative cost
  # false postive rate * number of negative instances * false positive cost +
  # false negative rate * number of positive instances * false negative cost
  cost <- perf@x.values[[1]] * cost.fp * sum( data[[actual]] == 0 ) +
    ( 1 - perf@y.values[[1]] ) * cost.fn * sum( data[[actual]] == 1 )
  
  cost_dt <- data.frame( cutoff = pred@cutoffs[[1]], cost = cost )
  
  # optimal cutoff value, and the corresponding true positive and false positive rate
  best_index  <- which.min(cost)
  best_cost   <- cost_dt[ best_index, "cost" ]
  best_tpr    <- roc_dt[ best_index, "tpr" ]
  best_fpr    <- roc_dt[ best_index, "fpr" ]
  best_cutoff <- pred@cutoffs[[1]][ best_index ]
  
  if (is.infinite(best_cutoff)) {
    best_cutoff = 0.5
  }
  
  f1 <- get_f1(data, predict, actual, best_cutoff)
  
  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # normalize the cost to assign colors to 1
  normalize <- function(v) ( v - min(v) ) / diff( range(v) )
  # create color from a palette to assign to the 100 generated threshold between 0 ~ 1
  # then normalize each cost and assign colors to it, the higher the blacker
  # don't times it by 100, there will be 0 in the vector
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)
  col_by_cost <- col_ramp[ ceiling( normalize(cost) * 99 ) + 1 ]
  
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) +
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.2 ) +
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) +
    labs( title = "ROC", x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = best_tpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = best_fpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" )
  
  cost_plot <- ggplot( cost_dt, aes( cutoff, cost ) ) +
    geom_line( color = "blue", alpha = 0.5 ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.5 ) +
    ggtitle( "Cost" ) +
    #  scale_y_continuous( labels = comma ) +
    geom_vline( xintercept = best_cutoff, alpha = 0.8, linetype = "dashed", color = "steelblue4" )
  
  options(scipen = '999')
  # the main title for the two arranged plot
  sub_title <- sprintf(other_title,  "Cutoff at %.2f - Total Cost = %a, AUC = %.3f",
                       best_cutoff, best_cost, auc )
  
  # arranged into a side by side plot
  plot <- arrangeGrob( roc_plot, cost_plot, ncol = 2,
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ) ) )
  
  # Calculate NPV
  cm <- table(data[[actual]], ifelse(data[[predict]] >= best_cutoff,1,0))
  npv <- cm[1,1]/(cm[1,1] + cm[2,1])
  
  # Calculate confidence intervals
  n <- nrow(data)
  boot_n <- 2000
  
  # Bootstrap for AUC CI with error handling
  auc_boot <- replicate(boot_n, {
    idx <- sample(1:n, n, replace=TRUE)
    # Check if we have both classes present
    if (length(unique(data[[actual]][idx])) == 2) {
      pred_boot <- prediction(data[[predict]][idx], data[[actual]][idx])
      tryCatch({
        performance(pred_boot, "auc")@y.values[[1]]
      }, error = function(e) NA)
    } else {
      NA  # Return NA if we don't have both classes
    }
  })
  # Remove any NA values before calculating quantiles
  auc_boot <- auc_boot[!is.na(auc_boot)]
  if (length(auc_boot) > 0) {
    auc_ci <- quantile(auc_boot, c(0.025, 0.975))
  } else {
    auc_ci <- c(NA, NA)
  }
  
  # Bootstrap for F1 CI with error handling
  f1_boot <- replicate(boot_n, {
    idx <- sample(1:n, n, replace=TRUE)
    if (length(unique(data[[actual]][idx])) == 2) {
      tryCatch({
        get_f1(data[idx,], predict, actual, best_cutoff)
      }, error = function(e) NA)
    } else {
      NA
    }
  })
  f1_boot <- f1_boot[!is.na(f1_boot)]
  f1_ci <- quantile(f1_boot, c(0.025, 0.975))
  
  # Bootstrap for NPV CI with error handling
  npv_boot <- replicate(boot_n, {
    idx <- sample(1:n, n, replace=TRUE)
    if (length(unique(data[[actual]][idx])) == 2) {
      cm_boot <- table(data[[actual]][idx], ifelse(data[[predict]][idx] >= best_cutoff,1,0))
      if (ncol(cm_boot) == 2 && nrow(cm_boot) == 2) {
        cm_boot[1,1]/(cm_boot[1,1] + cm_boot[2,1])
      } else {
        NA
      }
    } else {
      NA
    }
  })
  npv_boot <- npv_boot[!is.na(npv_boot)]
  npv_ci <- quantile(npv_boot, c(0.025, 0.975))
  
  return(list(data = data,
              plot = plot,
              pred = pred,
              perf = perf,
              cutoff = best_cutoff,
              auc = auc,
              sensitivity = best_tpr,
              specificity = 1 - best_fpr,
              f1 = f1,
              totalcost = best_cost,
              npv = npv,
              auc_ci = auc_ci,
              f1_ci = f1_ci,
              npv_ci = npv_ci))
}  

ROCInfo_atcutoff_ci <- function(data, predict, actual, cutoff, other_title, model) {
  # Handle model-specific prediction column naming
  if (model %in% c("enet", "rf","gbm","svmLinear","svmRadial","svmPoly","xgboost","nnet")) {
    predict <- paste0(predict,".Yes")
  }
  
  # Convert factors to numeric
  if (any(unique(data[[actual]]) %in% c("Yes","No"))) {
    data[[actual]] <- as.integer(data[[actual]] == "Yes")
  }
  data[[predict]] <- as.numeric(data[[predict]])
  
  # Calculate base confusion matrix with error handling
  pred_classes <- ifelse(data[[predict]] >= cutoff, 1, 0)
  cm <- table(factor(data[[actual]], levels=c(0,1)), 
              factor(pred_classes, levels=c(0,1)))
  
  # Ensure confusion matrix has proper dimensions
  if (nrow(cm) != 2 || ncol(cm) != 2) {
    cm <- matrix(0, nrow=2, ncol=2)
    dimnames(cm) <- list(c("0","1"), c("0","1"))
    # Fill available values
    actual_classes <- as.character(data[[actual]])
    pred_classes <- as.character(pred_classes)
    for(i in unique(actual_classes)) {
      for(j in unique(pred_classes)) {
        if(i %in% rownames(cm) && j %in% colnames(cm)) {
          cm[i,j] <- sum(actual_classes == i & pred_classes == j)
        }
      }
    }
  }
  
  # Extract confusion matrix values with safety checks
  tp <- cm["1","1"]
  tn <- cm["0","0"]
  fp <- cm["0","1"]
  fn <- cm["1","0"]
  
  # Calculate base metrics with safety checks
  sensitivity <- ifelse(tp + fn > 0, tp/(tp + fn), 0)
  specificity <- ifelse(tn + fp > 0, tn/(tn + fp), 0)
  precision <- ifelse(tp + fp > 0, tp/(tp + fp), 0)
  npv <- ifelse(tn + fn > 0, tn/(tn + fn), 0)
  f1 <- ifelse(precision + sensitivity > 0, 
               2 * (precision * sensitivity)/(precision + sensitivity), 
               0)
  
  # Calculate ROC and AUC
  pred <- prediction(data[[predict]], data[[actual]])
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  roc_dt <- data.frame(fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
  
  # Bootstrap confidence intervals
  n <- nrow(data)
  n_boot <- 2000
  boot_stats <- replicate(n_boot, {
    # Sample with replacement
    boot_idx <- sample(n, replace = TRUE)
    boot_data <- data[boot_idx,]
    
    # Calculate bootstrap prediction
    boot_pred <- prediction(boot_data[[predict]], boot_data[[actual]])
    
    # Calculate bootstrap metrics
    boot_cm <- table(
      Actual = boot_data[[actual]], 
      Predicted = ifelse(boot_data[[predict]] >= cutoff, 1, 0)
    )
    
    # Handle cases where confusion matrix doesn't have all cells
    if (nrow(boot_cm) < 2 || ncol(boot_cm) < 2) {
      return(c(sens = NA, spec = NA, f1 = NA, npv = NA, auc = NA))
    }
    
    # Calculate bootstrap values
    boot_tp <- boot_cm[2,2]
    boot_tn <- boot_cm[1,1]
    boot_fp <- boot_cm[1,2]
    boot_fn <- boot_cm[2,1]
    
    boot_sens <- boot_tp / (boot_tp + boot_fn)
    boot_spec <- boot_tn / (boot_tn + boot_fp)
    boot_prec <- boot_tp / (boot_tp + boot_fp)
    boot_npv <- boot_tn / (boot_tn + boot_fn)
    
    # Handle edge cases for F1
    if (boot_tp == 0) {
      boot_f1 <- 0
    } else {
      boot_f1 <- 2 * (boot_prec * boot_sens) / (boot_prec + boot_sens)
    }
    
    # Calculate AUC
    tryCatch({
      boot_auc <- performance(boot_pred, "auc")@y.values[[1]]
    }, error = function(e) {
      boot_auc <- NA
    })
    
    c(sens = boot_sens, 
      spec = boot_spec, 
      f1 = boot_f1,
      npv = boot_npv,
      auc = boot_auc)
  })
  
  # Initialize ci_ranges with default values
  ci_ranges <- matrix(NA, nrow=5, ncol=2,
                     dimnames=list(c("sens", "spec", "f1", "npv", "auc"),
                                 c("2.5%", "97.5%")))
  
  # Safely calculate confidence intervals
  if(!is.null(dim(boot_stats))) {
    for(metric in rownames(ci_ranges)) {
      metric_values <- boot_stats[metric,]
      metric_values <- metric_values[!is.na(metric_values)]
      if(length(metric_values) > 0) {
        ci_ranges[metric,] <- quantile(metric_values, probs=c(0.025, 0.975))
      }
    }
  }

  # Create plot
  roc_plot <- ggplot(roc_dt, aes(fpr, tpr)) + 
    geom_line(color = rgb(0, 0, 1, alpha = 0.3)) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), alpha = 0.8, color = "royalblue") + 
    labs(title = "ROC", x = "False Positive Rate", y = "True Positive Rate") +
    geom_hline(yintercept = sensitivity, alpha = 0.8, linetype = "dashed", color = "steelblue4") +
    geom_vline(xintercept = 1-specificity, alpha = 0.8, linetype = "dashed", color = "steelblue4")
  
  # Create plot title with proper formatting
  options(scipen = '999')
  sub_title <- paste("Cutoff =", round(cutoff, 3), "AUC =", round(auc, 3))
  plot <- arrangeGrob(roc_plot, ncol = 2, 
                      top = textGrob(sub_title, gp = gpar(fontsize = 16, fontface = "bold")))
  
  return(list(
    data = data,
    plot = plot,
    pred = pred,
    perf = perf,
    cutoff = cutoff,
    auc = auc,
    sensitivity = sensitivity,
    specificity = specificity,
    f1 = f1,
    npv = npv,
    cm = cm,
    sens_ci = ci_ranges["sens",],
    spec_ci = ci_ranges["spec",],
    f1_ci = ci_ranges["f1",],
    npv_ci = ci_ranges["npv",],
    auc_ci = ci_ranges["auc",]
  ))
}

ConfusionMatrixInfo <- function( data, predict, actual, cutoff, get_plot, data_type, points,model) {	
  if (model  %in% c("rf","svmLinear","gbm","xgboost","svmRadial","enet","nnet","svmPoly") ) {
    predict <- paste0(predict,".Yes")
  } 
  
  if (any(!unique(data[[actual]]) %in% c("Yes","No"))) {
    data[[actual]] <- ifelse(data[[actual]] == 1,"Yes","No")
  }
  
  # extract the column ;
  # relevel making 1 appears on the more commonly seen position in 
  # a two by two confusion matrix	
  predict <- data[[predict]]
  actual  <- factor( data[[actual]], levels=c("Yes","No"))
  
  if(data_type == 'null') {
    
    age <- data$agesamplecollection
    show_legend <- TRUE
    legend_cols <- c('blue', 'purple')  # FP in blue, TN in purple
    breaks <- c("FP", "TN")
    
    
  } else {
    age <- data$ageofonset
    show_legend <- TRUE
    legend_cols <- c('red', 'orange', 'blue', 'purple')  # TN in purple, FP in blue
    breaks <- c( "TP", "FN", "FP", "TN" )
  }
  result <- data.table( actual = actual, predict = predict, age = age)
  
  # caculating each pred falls into which category for the confusion matrix
  result[ , type := ifelse( predict >= cutoff & actual == 'Yes', "TP",
                            ifelse( predict >= cutoff & actual == 'No', "FP", 
                                    ifelse( predict <  cutoff & actual == 'Yes', "FN", "TN" ))) %>% 
            as.factor() ]
  
  result <- as.data.frame(result)
  if(data_type == 'null'){
    result$actual <- as.factor(ifelse(result$actual == 'Yes', 'Cancer before 6', 'Cancer after 6'))
    result$acual <- factor(result$actual, levels = c('Cancer before 6', 'Cancer after 6'))
    x_lab = 'Age of Sample Collection'
  } else {
    result$actual <- as.factor(ifelse(result$actual == 'Yes', 'Cancer before 6', 'Cancer after 6'))
    result$acual <- factor(result$actual, levels = c('Cancer before 6', 'Cancer after 6'))
    x_lab = 'True Age of Cancer Onset'
    
  }
  result$actual<- with(result, reorder(actual, acual, function(x) length(x)))  
  library(ggthemes)
  library(ggrepel)
  # jittering : can spread the points along the x axis 
  result$age <- round(result$age/12)
  if(points){
    plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
      geom_jitter(size = 3, alpha = 0.8, show.legend = show_legend) +
      scale_color_manual(name = '',
                         values = legend_cols,
                         breaks = breaks)+
      geom_violin( fill = "darkgrey", color = NA , alpha = 0.3) +
      geom_hline( yintercept = cutoff, color = 'black', alpha = 0.6, linetype = 2 ) + 
      geom_vline(xintercept = 1.5, linetype = 2) +
      scale_y_continuous( limits = c( 0, 1 ) ) + 
      geom_text(aes(label = ifelse(type %in% c("FP","FN"),age,'')),alpha = 0.7, cex = 5,position=position_jitter(width = 0.4, height = 0), show.legend = FALSE)+
      #      ggtitle( sprintf( other_title,"_Cutoff at %.2f", cutoff ) ) +
      theme(legend.position="bottom", text = element_text(size=16)) +
      xlab(x_lab) + ylab('Predicted Probability of Early Cancer Onset') + theme_bw(base_size = 18) +
      theme(panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
  } else {
    plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
      geom_point(size = 0, show.legend = FALSE) +
      scale_color_manual(name = 'Result',
                         values = c('red', 'orange','green', 'blue'),
                         breaks = c( "TP", "FN", "FP", "TN" ))+
      geom_violin( fill = "lightgrey", color = NA) +
      geom_hline( yintercept = cutoff, color = 'black', alpha = 0.6, linetype = 2 ) + 
      geom_text(aes(label = age),alpha = 0.7, cex = 5,position=position_jitter(width = 0.4, height = 0), show.legend = TRUE)+
      geom_vline(xintercept = 1.5, linetype = 2) +
      scale_y_continuous( limits = c( 0, 1 ) ) + 
      guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
      #    ggtitle( sprintf( other_title,"_Cutoff at %.2f", cutoff ) ) +
      theme(text = element_text(size=14)) + 
      xlab(x_lab) + ylab('Predicted Probability of Early Cancer Onset')  + theme_bw(base_size = 18) +
      theme(panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
  }
  
  
  
  if(get_plot) {
    return(plot)
  } else {
    return(as.data.frame(result))
  }
}

reliability.plot <- function(obs, pred, bins=10, scale=T) {
  # Plots a reliability chart and histogram of a set of predicitons from a classifier
  #
  # Args:
  # obs: Vector of true labels. Should be binary (0 or 1)
  # pred: Vector of predictions of each observation from the classifier. Should be real
  # number
  # bins: The number of bins to use in the reliability plot
  # scale: Scale the pred to be between 0 and 1 before creating reliability plot
  require(plyr)
  library(Hmisc)
  min.pred <- min(pred)
  max.pred <- max(pred)
  min.max.diff <- max.pred - min.pred
  if (scale) {
    pred <- (pred - min.pred) / min.max.diff
  }
  bin.pred <- cut(pred, bins)
  k <- ldply(levels(bin.pred), function(x) {
    idx <- x == bin.pred
    c(sum(obs[idx]) / length(obs[idx]), mean(pred[idx]))
  })
  is.nan.idx <- !is.nan(k$V2)
  k <- k[is.nan.idx,]
  return(k)
}

platt_scaling <- function(train_results, test_results,model) {
  actual <- "test_label"
  if (model %in% c("rf","svmLinear","xgboost","svmRadial","gbm","nnet") ) {
    predict <- "test_pred.Yes"
    #    tmp <- ifelse(train_results[[actual]] == "Yes", 1, 0)
  } else {
    predict <- "test_pred"
  }
  #calculating Log Loss without Platt Scaling
  ll <- LogLoss(train_results[[predict]],train_results[[actual]])
  ll_df <- data.frame(x=train_results[[predict]],y=as.factor(train_results[[actual]]))
  model_log <-glm(y~x,data = ll_df,family = binomial)
  #predicting on the cross validation after platt scaling
  result_platt<-predict(model_log,ll_df["x"],type = "response")
  train_results$test_pred_calibrated.Yes <- result_platt
  
  plot(c(0,1),c(0,1), col="grey",type="l",xlab = "Mean Prediction",ylab="Observed Fraction")
  # The line below computes the reliability plot data for cross validation dataset without platt scaling
  k <-reliability.plot(train_results[[actual]], train_results[[predict]],bins = 5)
  lines(k$V2, k$V1, xlim=c(0,1), ylim=c(0,1), xlab="Mean Prediction", ylab="Observed Fraction", col="red", type="o", main="Reliability Plot")
  
  #This line below computes the reliability plot data for cross validation dataset with platt scaling
  k <-reliability.plot(train_results[[actual]],result_platt,bins = 5)
  lines(k$V2, k$V1, xlim=c(0,1), ylim=c(0,1), xlab="Mean Prediction", ylab="Observed Fraction", col="blue", type="o", main="Reliability Plot")
  
  legend("topright",lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"),legend = c("platt scaling","without platt scaling"))
  
  # Predicting on the test dataset using Platt Scaling
  ll_df_test<-data.frame(x=test_results[[predict]])
  result_test_platt<-predict(model_log,ll_df_test,type="response")
  test_results$test_pred_calibrated.Yes <- result_test_platt
  
  plot(c(0,1),c(0,1), col="grey",type="l",xlab = "Mean Prediction",ylab="Observed Fraction")
  # The line below computes the reliability plot data for cross validation dataset without platt scaling
  k <-reliability.plot(test_results[[actual]], test_results[[predict]],bins = 5)
  lines(k$V2, k$V1, xlim=c(0,1), ylim=c(0,1), xlab="Mean Prediction", ylab="Observed Fraction", col="red", type="o", main="Reliability Plot")
  
  #This line below computes the reliability plot data for cross validation dataset with platt scaling
  k <-reliability.plot(test_results[[actual]],result_test_platt,bins = 5)
  lines(k$V2, k$V1, xlim=c(0,1), ylim=c(0,1), xlab="Mean Prediction", ylab="Observed Fraction", col="blue", type="o", main="Reliability Plot")
  
  legend("topright",lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"),legend = c("platt scaling","without platt scaling"))
  
  return(list(train_results,test_results))
}

plot.scores.AUC <- function(y, y.hat, measure = "tpr", x.measure = "fpr") {
  # y <- ifelse(y == "Yes",1,0)
  y <- as.numeric(y)
  par(mfrow=c(1,2))
  hist(y.hat[y == 1], col=rgb(1,0,0,0.5), 
       main = "Score Distribution",
       breaks=seq(min(y.hat),max(y.hat)+0.1, 0.05), xlab = "Prediction")
  hist(y.hat[y == 2], col = rgb(0,0,1,0.5), add=T, 
       breaks=seq(min(y.hat),max(y.hat)+0.1 , 0.05))
  legend("topright", legend = c("Class 0", "Class 1"),  col=c("red", "blue"), lty=1, cex=0.6)
  # plot ROC curve
  pr <- prediction(y.hat, y)
  prf <- performance(pr, measure = measure, x.measure = x.measure)
  # get AUC
  auc <- performance(pr, measure = "auc")@y.values[[1]]
  plot(prf, main = paste0("Curve (AUC: ", round(auc, 2), ")"))
}

dir <- "/Users/vallijahsubasri/Documents/lfs_ageofonset/data/rds/"
dir <- "/Users/vallijahsubasri/Documents/lfs_ageofonset/predictions/"

bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_72_FNW1_svmRadial"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_tp53_72_FNW1_svmRadial"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_tp53_clinicalonly_72_FNW1_rf"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_systreat_72_FNW1_xgboost"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_sex_tp53_72_PCA_FNW1_svmRadial"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_sex_tp53_72_FNW1_svmRadial"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_cellprops_sex_tp53_72_nometh_FNW1_svmRadial"
tp53_only <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_tp53_clinicalonly_72_FNW1_rf"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_family_tp53_72_FNW1_svmRadial"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_tp53_72_FNW1_svmRadial"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_sex_tp53_PCA_72_FNW1_svmRadial"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_sex_tp53_PCA_72_FNW1_xgboost"
bestmodel_name <- "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_cellprops_72_FNW1_svmRadial"

model_type = "svmRadial"
model_name <- unlist(str_split(bestmodel_name,"_"))[length(unlist(str_split(bestmodel_name,"_")))]

ROCobj_val <- readRDS(paste0(dir,bestmodel_name,'_ageofonset_ROCInfoVal.rds'))
val_results <- ROCobj_val[[1]]
if ("Yes" %in% val_results$test_label) {
  val_results$test_label <- ifelse(val_results$test_label=="Yes",1,0)
}
val_results <- val_results[!duplicated(val_results$ids),]
#val_results <- val_results[!val_results$tm_donor %in% c("4475"),]
#val_results <- val_results[!is.na(val_results$agesamplecollection),]
val_pr<-pr.curve(scores.class0 = val_results$test_pred.Yes[val_results$test_label == 0], scores.class1 = val_results$test_pred.Yes[val_results$test_label == 1])
ROCobj_val[[10]] <- val_pr$auc.integral
weight=2
pred_col='test_pred_calibrated'
ROCobj_val <- ROCInfo_atopt_ci(val_results,pred_col,"test_label",1,weight,bestmodel_name, model_type)
cutoff <- ROCobj_val[[5]]
# For validation results
cat(paste0("Validation Set Metrics:\n"))
cat(paste0("Cutoff: ", round(ROCobj_val[[5]], 3), "\n"))
cat(paste0("AUC: ", round(ROCobj_val[[6]], 3), 
           " (95% CI: ", round(ROCobj_val$auc_ci[1], 3), "-", round(ROCobj_val$auc_ci[2], 2), ")\n"))
cat(paste0("Sensitivity: ", round(ROCobj_val[[7]], 3), "\n"))
cat(paste0("Specificity: ", round(ROCobj_val[[8]], 3), "\n"))
cat(paste0("F1 Score: ", round(ROCobj_val[[9]], 3),
           " (95% CI: ", round(ROCobj_val$f1_ci[1], 3), "-", round(ROCobj_val$f1_ci[2], 3), ")\n"))
cat(paste0("NPV: ", round(ROCobj_val$npv, 3),
           " (95% CI: ", round(ROCobj_val$npv_ci[1], 3), "-", round(ROCobj_val$npv_ci[2], 3), ")\n"))
cat(paste0("AUPRC: ", round(ROCobj_val[[10]], 3), "\n\n"))


points=TRUE
val_cancer <- ConfusionMatrixInfo(val_results[val_results$cancerstatus!="Unaffected",],pred_col,'test_label',cutoff,TRUE,'not null',points,model_name)
val_null <- ConfusionMatrixInfo(val_results[val_results$cancerstatus=="Unaffected",],pred_col,'test_label',cutoff,TRUE,'null',points,model_name)
val_cancer ; val_null

ROCobj_test <- readRDS(paste0(dir,bestmodel_name,'_ageofonset_ROCInfoTest.rds')) 
test_results <- ROCobj_test[[1]]
test_results <- test_results[!duplicated(test_results$ids),]
ROCobj_test <- ROCInfo_atcutoff_ci(test_results,pred_col,"test_label",cutoff,bestmodel_name,model_type)
test_pr<-pr.curve(scores.class0 = test_results$test_pred_calibrated.Yes[test_results$test_label == 0], scores.class1 = test_results$test_pred_calibrated.Yes[test_results$test_label == 1])
ROCobj_test[[10]] <- test_pr$auc.integral

# For test results 
cat(paste0("Test Set Metrics:\n"))
cat(paste0("AUC: ", round(ROCobj_test[[6]], 3),
           " (95% CI: ", round(ROCobj_test$auc_ci[1], 3), "-", round(ROCobj_test$auc_ci[2], 2), ")\n"))
cat(paste0("Sensitivity: ", round(ROCobj_test[[7]], 3), "\n"))
cat(paste0("Specificity: ", round(ROCobj_test[[8]], 3), "\n"))
cat(paste0("F1 Score: ", round(ROCobj_test[[9]], 3),
           " (95% CI: ", round(ROCobj_test$f1_ci[1], 3), "-", round(ROCobj_test$f1_ci[2], 3), ")\n"))
cat(paste0("NPV: ", round(ROCobj_test$npv, 3),
           " (95% CI: ", round(ROCobj_test$npv_ci[1], 3), "-", round(ROCobj_test$npv_ci[2], 3), ")\n"))




test_cancer <- ConfusionMatrixInfo(test_results[test_results$cancerstatus!="Unaffected",],pred_col,'test_label',cutoff,TRUE,'not null',points,model_name)
test_null <- ConfusionMatrixInfo(test_results[test_results$cancerstatus=="Unaffected",],pred_col,'test_label',cutoff,TRUE,'null',points,model_name)
test_cancer ; test_null

## Plot ROC curve
par(pty = "s")
pROC_val <- roc(val_results$test_label, val_results$test_pred_calibrated.Yes,
                # arguments for ci
                ci=TRUE, ci.alpha=0.95, stratified=FALSE,
                # arguments for plot
                plot=TRUE, grid=TRUE, show.thres=TRUE)
sens.ci <- ci.se(pROC_val)
plot(sens.ci, type="shape", col="lightgrey",axes=FALSE)
plot(sens.ci, type="bars", add=FALSE)
best.coords <- coords(pROC_val, "best", best.method="youden")
best.coords[1]=cutoff
abline(v=best.coords[1,'threshold'], lty=2, col="red")
#abline(h=best.coords["specificity"], lty=2, col="red")
abline(h=best.coords[1,"sensitivity"], lty=2, col="blue")
text(0.2, 0.4, paste("AUC:", round(pROC_val$auc, 3)), cex=0.75)
text(0.2, 0.35, paste("95% CI:", round(pROC_val$ci, 3)[1],"-" ,round(pROC_val$ci, 3)[3]), cex=0.75)

## Plot ROC curve
par(pty = "s")
pROC_test <- roc(test_results$test_label, test_results$test_pred_calibrated.Yes,
                 # arguments for ci
                 ci=TRUE, ci.alpha=0.95, stratified=FALSE,
                 # arguments for plot
                 plot=TRUE, grid=TRUE, show.thres=TRUE)
sens.ci <- ci.se(pROC_test,boot.n=10)
plot(sens.ci, type="shape", col="lightgrey",axes=FALSE)
plot(sens.ci, type="bars", add=FALSE)
best.coords <- coords(pROC_val, "best", best.method="youden")
best.coords[1]=cutoff
abline(v=best.coords[1,"threshold"], lty=2, col="red")
#abline(h=best.coords["specificity"], lty=2, col="red")
abline(h=best.coords[1,"sensitivity"], lty=2, col="blue")
text(0.2, 0.4, paste("AUC:", round(pROC_test$auc, 3)), cex=0.75)
text(0.2, 0.35, paste("95% CI:", round(pROC_test$ci, 3)[1],"-" ,round(pROC_test$ci, 3)[3]), cex=0.75)

## plot immune proportions

data <- readRDS("/Users/vallijahsubasri/Documents/lfs_ageofonset/data/NoobCorrected_beta_ProjPC2Adj_lfs_3UTR.rds")
#fwrite(as.list(colnames(data)[68:length(data)]),'~/Downloads/genes',sep='\n',quote=FALSE)
u <- umap(data[68:length(data)],n_neighbors=20)
#u <- umap(data[c("CLDN10", "CLDN18", "CTNNA1", "CTNNA3", "PIK3R1", "PXN", "RASSF5")])
u <- cbind(data.frame(u$layout),data[1:67])
u$earlyonset <- ifelse((u$ageofonset > 72 | is.na(u$ageofonset)),"No","Yes")
ggplot(u, aes(X1,X2,color=ageofsamplecollection)) +
  geom_point() + 
  theme_minimal()

n <- length(names(table(u$cancer_diagnosis))[table(u$cancer_diagnosis)>4])
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(u[u$cancer_diagnosis %in% names(table(u$cancer_diagnosis))[table(u$cancer_diagnosis)>4],], aes(X1,X2,color=cancer_diagnosis)) +
  geom_point() + 
  theme_minimal() + 
  scale_color_manual(values=col_vector)

n <- length(names(table(u$tissue_type))[table(u$tissue_type)>1])
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

tissue_types=names(table(u$tissue_type))[table(u$tissue_type)>1]
ggplot(u[u$tissue_type %in% tissue_types,], aes(X1,X2,color=tissue_type)) +
  geom_point(size=2) + 
  theme_minimal() + 
  scale_color_manual(values=col_vector)

p <- prcomp(data[68:length(data)])
p <- cbind(data.frame(p$x),data[1:67])
p$earlyonset <- ifelse((p$ageofonset > 72 | is.na(p$ageofonset)),"No","Yes")

ggplot(p, aes(PC1,agesamplecollection,color=agesamplecollection)) +
  geom_point() + 
  theme_minimal() +
  stat_cor(label.y.npc="top", label.x.npc = "left", method = "pearson",size=2.5)+
  labs(y= "Age of Sample Collection (months) ", x = "PC1")  + 
  scale_color_continuous(name = "Age of Sample Collection")

TP <- c(val_results$ids[val_results$test_label == 1 & val_results$test_pred_calibrated.Yes > cutoff],
        test_results$ids[test_results$test_label == 1 & test_results$test_pred_calibrated.Yes > cutoff])
FN <- c(val_results$ids[val_results$test_label == 1 & val_results$test_pred_calibrated.Yes < cutoff],
        test_results$ids[test_results$test_label == 1 & test_results$test_pred_calibrated.Yes < cutoff])


test_positive <- test_results[test_results$test_label == 1,] ; test_positive$label <- ifelse(test_positive$test_pred_calibrated.Yes > cutoff,"TP","FN")
val_pos_agediff <- ggplot(test_positive,aes(test_pred_calibrated.Yes,age_diff)) + 
  geom_point(aes(color=label),size=2) + 
  geom_smooth(color="light grey",alpha=0.2,method='lm', se=TRUE, fullrange=FALSE, level=0.95)+
  stat_cor(method = "pearson",label.x.npc="middle", label.y.npc="top",family="Helvetica Neue")+
  theme_classic() + 
  theme(text = element_text(size=12,family="Helvetica Neue")) +
  labs(x ="True Class Predicted Probability", y="|Age at Sample collection - Age at Diagnosis|") + 
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  guides(color=guide_legend(title="Predicted Outcome"))

val_positive <- val_results[val_results$test_label == 1,] ; val_positive$label <- ifelse(val_positive$test_pred_calibrated.Yes > cutoff,"TP","FN")
test_pos_agediff <- ggplot(val_positive,aes(test_pred_calibrated.Yes,age_diff)) + 
  geom_point(aes(color=label),size=2) + 
  geom_smooth(color="light grey",alpha=0.2,method='lm', se=TRUE, fullrange=FALSE, level=0.95)+
  stat_cor(method = "pearson",label.x.npc="middle", label.y.npc="top",family="Helvetica Neue")+
  theme_classic() +
  theme(text = element_text(size=12,family="Helvetica Neue")) +
  labs(x ="True Class Predicted Probability", y="|Age at Sample collection - Age at Diagnosis|") + 
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  guides(color=guide_legend(title="Predicted Outcome"))


ggarrange(val_pos_agediff,
          test_pos_agediff,
          nrow=2,
          labels=c("A","B"),
          common.legend = TRUE,
          legend="bottom")
## plot clinical variables
train_dat <- readRDS(paste0(dir,'NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_TrainingSet.rds')) 
train_dat$dataset <- "train"
test_dat <- readRDS(paste0(dir,'NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_TestSet.rds')) 
test_dat$dataset <- "val"
val_dat <- readRDS(paste0(dir,'NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_ValidationSet.rds')) 
val_dat$dataset <- "test"
alldata <- rbind(train_dat,test_dat,val_dat)

ggplot(alldata[!is.na(alldata$cancer_atdraw),],aes(x=cancer_atdraw,..count..,fill=cancerstatus)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~dataset) + 
  theme_classic() 

ggplot(alldata[!is.na(alldata$cancer_atdraw),],aes(x=cancer_atdraw,..count..,fill=cancerstatus)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~dataset) + 
  theme_classic()

