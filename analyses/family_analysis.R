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

ROCInfo_atopt <- function( data, predict, actual, cost.fp, cost.fn, other_title, model) {
  
  if (model %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet") ) {
    predict <- paste0(predict,".Yes")
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
  
  return( list(data = data,
               plot      = plot,
               pred      = pred,
               perf      = perf,
               cutoff    = best_cutoff,
               auc         = auc,
               sensitivity = best_tpr,
               specificity = 1 - best_fpr,
               f1 = f1,
               totalcost   = best_cost) )
}  

ROCInfo_atcutoff <- function( data, predict, actual, cutoff, other_title, model) {
  
  if (model %in% c("enet","rf","gbm","svmLinear","svmRadial","xgboost","nnet") ) {
    predict <- paste0(predict,".Yes")
  }
  
  if (any(unique(data[[actual]]) %in% c("Yes","No"))) {
    data[[actual]] <- ifelse(data[[actual]] == "Yes", 1, 0)
  }
  
  # calculate the values using the ROCR library
  # true positive, false postive 
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  ## get sens/spec/f1
  cm <- table(data[[actual]], ifelse(data[[predict]] >= cutoff,1,0) )
  if (!("1" %in% colnames(cm))) {
    cm <- cbind(cm,c(0,0))
  } else if (!("0" %in% colnames(cm))) {
    cm <- cbind(c(0,0),cm)
  }
  sensitivity <- cm[2,2]/sum(cm[2,])
  specificity <- cm[1,1]/sum(cm[1,])
  accuracy <- (cm[1,1] + cm[2,2])/sum(cm)
  f1 <- get_f1(data, predict, actual, cutoff)
  
  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # create color from a palette to assign to the 100 generated threshold between 0 ~ 1
  # then normalize each cost and assign colors to it, the higher the blacker
  # don't times it by 100, there will be 0 in the vector
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)
  
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) +
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) +
    labs( title = "ROC", x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = sensitivity, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = 1-specificity, alpha = 0.8, linetype = "dashed", color = "steelblue4" )
  
  options(scipen = '999')
  # the main title for the two arranged plot
  sub_title <- sprintf(other_title,  "Cutoff at %.2f , AUC = %.3f",
                       cutoff, auc )
  
  # arranged into a side by side plot
  plot <- arrangeGrob( roc_plot, ncol = 2,
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ) ) )
  
  return( list( data = data,
                plot      = plot,
                pred      = pred,
                perf      = perf,
                cutoff    = cutoff,
                auc         = auc,
                sensitivity = sensitivity,
                specificity = specificity,
                f1 = f1) )
}

ConfusionMatrixInfo <- function( data, predict, actual, cutoff, get_plot, data_type, points,model) {	
  if (model  %in% c("rf","svmLinear","gbm","xgboost","svmRadial","enet","nnet") ) {
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
    legend_cols <- c('blue', 'blue', 'blue', 'blue')
    breaks <- c("TN" )
    
    
  } else {
    age <- data$ageofonset
    show_legend <- TRUE
    legend_cols <- c('red', 'orange','blue', 'purple')
    breaks <- c( "TP", "FN", "FP", "TN" )
  }
  result <- data.table( actual = actual, predict = predict, age = age)
  
  # caculating each pred falls into which category for the confusion matrix
  result[ , type := ifelse( predict >= cutoff & actual == 'Yes', "TP",
                            ifelse( predict >= cutoff & actual == 'No', "FP", 
                                    ifelse( predict <  cutoff & actual == 'Yes', "FN", "TN" ) ) ) %>% as.factor() ]
  
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
  library(ROCR)
  pr <- prediction(y.hat, y)
  prf <- performance(pr, measure = measure, x.measure = x.measure)
  # get AUC
  auc <- performance(pr, measure = "auc")@y.values[[1]]
  plot(prf, main = paste0("Curve (AUC: ", round(auc, 2), ")"))
}

dir <- "/Users/vallijahsubasri/Documents/lfs_ageofonset/family/seed2/"

model_name <- "svmRadial"
ROCobj_val <- readRDS(paste0(dir,"NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_family_sex_72_FNW1_svmRadial_ageofonset_ROCInfoTest.rds"))
ROCobj_test <- readRDS(paste0(dir,"NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_family_sex_72_FNW1_svmRadial_ageofonset_ROCInfoVal.rds"))
model <- readRDS(paste0(dir, "NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_family_sex_72_svmRadial_ageofonset_results.rds"))

val_results <- ROCobj_val[[1]]
val_results <- val_results[!duplicated(val_results$ids),]
val_results <- val_results[!val_results$tm_donor %in% c("4475"),]
val_pr<-pr.curve(scores.class0 = val_results$test_pred.Yes[val_results$test_label == 0], scores.class1 = val_results$test_pred.Yes[val_results$test_label == 1])
ROCobj_val[[10]] <- val_pr$auc.integral
weight=1
pred_col='test_pred_calibrated'
ROCobj_val <- ROCInfo_atopt(val_results,pred_col,"test_label",1,weight,bestmodel_name, model_type)
cutoff <- 0.4#ROCobj_val[[5]]
cat(paste0("Cutoff: ",round(ROCobj_val[[5]],2)))
cat(paste0("Auc: ",round(ROCobj_val[[6]],2)))
cat(paste0("Sensitivity: ",round(ROCobj_val[[7]],2)))
cat(paste0("Specificity: ",round(ROCobj_val[[8]],2)))
cat(paste0("F1 Score: ",round(ROCobj_val[[9]],3)))

points=TRUE
val_cancer <- ConfusionMatrixInfo(val_results[val_results$cancerstatus!="Unaffected",],pred_col,'test_label',cutoff,TRUE,'not null',points,model_name)
val_null <- ConfusionMatrixInfo(val_results[val_results$cancerstatus=="Unaffected",],pred_col,'test_label',cutoff,TRUE,'null',points,model_name)
val_cancer ; val_null

test_results <- ROCobj_test[[1]]
test_results <- test_results[!duplicated(test_results$ids),]
ROCobj_test <- ROCInfo_atcutoff(test_results,pred_col,"test_label",cutoff,bestmodel_name,model_type)
test_pr<-pr.curve(scores.class0 = test_results$test_pred_calibrated.Yes[test_results$test_label == 0], scores.class1 = test_results$test_pred_calibrated.Yes[test_results$test_label == 1])
ROCobj_test[[10]] <- test_pr$auc.integral
cat(paste0("Auc: ",round(ROCobj_test[[6]],2)))
cat(paste0("Sensitivity: ",round(ROCobj_test[[7]],2)))
cat(paste0("Specificity: ",round(ROCobj_test[[8]],2)))
cat(paste0("F1 Score: ",round(ROCobj_test[[9]],2)))

test_cancer <- ConfusionMatrixInfo(test_results[test_results$cancerstatus!="Unaffected",],pred_col,'test_label',cutoff,TRUE,'not null',points,model_name)
test_null <- ConfusionMatrixInfo(test_results[test_results$cancerstatus=="Unaffected",],pred_col,'test_label',cutoff,TRUE,'null',points,model_name)
test_cancer ; test_null

