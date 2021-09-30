#####################################################################
############################# Libraries #############################
#####################################################################

require(reshape2)
require(minfi)
require(dplyr)
require(doParallel)
require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
require(stringr)
require(ggplot2)
require(ROCR)
require(data.table)
require(caret)
require(gridExtra)
require(grid)
require(MLmetrics)

#####################################################################
############################# Functions #############################
#####################################################################

##########
# Function to remove batch effect
#'
#' @rdname remove_batch_confounder
#' @name remove_batch_confounder
#' 
#' @param data Dataframe containing methylation data and clinical variables prior to batch confounder removal
#' @return Dataframe containing methylation data and clinical variables after batch confounder removal
##########

remove_batch_confounder <- function(data) {
  p_train <- readRDS('rds/batchEffectRemoval.rds')
  cols <- grepl( "cg" , names( data ))
  beta <- data[ , cols ]
  clin <- data[ , -cols ]
  test_pred <-predict(p_train,beta)
  Xhat_pred <- test_pred[, !(colnames(test_pred) %in% c(maxpc))] %*% t(p_train$rotation[, !(colnames(p_train$rotation) %in% c(maxpc))])
  beta_adj_pred <- scale(Xhat_pred, center = -(colMeans(beta)), scale = T)
  data_adj_test <- cbind(clin,beta_adj_pred)
  return(data_adj_test)
}

##########
# Function to batch 450k data onto the 850k space
#'
#' @rdname remove_array_confounder
#' @name remove_array_confounder
#' 
#' @param data Dataframe containing methylation data and clinical variables prior to array confounder removal
#' @return Dataframe containing methylation data and clinical variables after array confounder removal
##########

remove_array_confounder <- function(data) {
  cat("[  450k data onto 850k space ]","\n")
  models <- readRDS('rds/arrayEffectRemoval.rds')
  probes <- grepl( "cg" , names( data )); lm_all <- list()
  for (probe in probes) {
    probe_model=models[[probe]]
    new_data <- data.frame(probe_450 = data[,probe],
         age_sample_collection = data$agesamplecollection,
         gender = data$gender)
    lm_probe <- predict(probe_model,new_data)
    if (i %% 50000 == 0) {
      print(i)
    }
  }
  lm_all[[probe]] <-lm_probe
  corrected_meth <- do.call(cbind,lm_all)
  corrected_data <- cbind(data[,-cols],corrected_meth)
  return(corrected_data)
}

##########
# Function that identifies ids of outliers
#'
#' @rdname get_outliers
#' @name get_outliers
#' 
#' @param pc Dataframe containing first two principal components, must contain columnns: PC1, PC2, SentrixID
#' @param n Number of standard deviations as lower and upper bound for outlier removal
#' @return List of ids to keep after outlier removal
##########

get_outliers <- function(pc,n) {
  pc$PC1 <- scale(pc$PC1)
  pc$PC2 <- scale(pc$PC2)
  u_thres_pc1 <- mean(pc$PC1) + n*sd(pc$PC1)
  l_thres_pc1 <- mean(pc$PC1) - n*sd(pc$PC1)
  u_thres_pc2 <- mean(pc$PC2) + n*sd(pc$PC2)
  l_thres_pc2 <- mean(pc$PC2) - n*sd(pc$PC2)
  pc2 <- pc[pc$PC1 < u_thres_pc1 & pc$PC1 > l_thres_pc1, ]
  pc2 <- pc2[pc2$PC2 < u_thres_pc2 & pc2$PC2 > l_thres_pc2, ]
  outlier_id <- as.character(pc$SentrixID[!pc$SentrixID %in% pc2$SentrixID])
  cat(paste0(length(outlier_id)," outliers removed","\n"))
  keep_id <- as.character(pc2$SentrixID)
  return(keep_id)
}

##########
# Function that quantifies association between PCs and confounders
#'
#' @rdname generate_pcsummary
#' @name generate_pcsummary
#' 
#' @param pc_clin Dataframe with principal components (i.e. PC1, PC2, PC3, etc.) and confounder variables: array, batch, agesamplecollection
#' @return Dataframe containing association between principal components and confounders
##########

generate_pcsummary <- function(pc_clin) {
  pcs <- colnames(pc_clin)[grepl( "PC" , names(pc_clin) )] ;
  allp_array <- list() ; allp_batch <- list() ; allp_agesamplecollection<- list()
  for (i in pcs) {
    allp_array[[i]] <- wilcox.test(pc_clin[,i] ~ as.factor(pc_clin$array))$p.value
    allp_batch[[i]] <- summary(aov(pc_clin[,i] ~ as.factor(pc_clin$batch)))[[1]][1,"Pr(>F)"]
    allp_batch[[i]] <- cor.test(pc_clin[,i] ~ as.factor(pc_clin$agesamplecollection)$p.value
  }

  allp <- data.frame(PC=pcs,
                      batch_p=as.numeric(do.call("cbind",allp_batch)),
                      array_p=as.numeric(do.call("cbind",allp_array)),
                      agesamplecollection_p=as.numeric(do.call("cbind",allp_agesamplecollection)))
  allp <- cbind(allp,data.frame(t(summary(pc)$importance)))
  allp$batch_padj <- allp$batch_p/allp$Proportion.of.Variance
  allp$array_padj <- allp$array_p/allp$Proportion.of.Variance
  allp$agesamplecollection_padj <- allp$agesamplecollection_p/allp$Proportion.of.Variance
  return(allp)
}


##########
# Function to scale methylation variables
#'
#' @rdname scale_df
#' @name scale_df
#' 
#' @param data Dataframe containing unscaled methylation data and clinical variables
#' @param genes List containing gene-wise methylation features to scale
#' @return Dataframe containing scaled methylation data and clinical variables
##########

scale_df <- function(data,genes) {
  tmp <- data[genes]
  tmp <- scale(tmp, center=TRUE,scale=TRUE)
  data_scaled <- cbind(data[,-genes],tmp)
  return(data_scaled)
}

##########
# Function to extract probes based on location and aggregate by gene
#'
#' @rdname aggregate_probes
#' @name aggregate_probes
#' 
#' @param data Dataframe containg all the methylation data and clinical variables
#' @param features Dataframe containing probe and gene-wise features
#' @return Dataframe containg all the methylation data and clinical variables
##########

aggregate_probes <- function(data,features) {
  allprobes <- colnames(data)[grepl('cg',colnames(data))]
  clincols <- colnames(data)[-allprobes]
  feat_data <- data[c(clincols,features$probe)]
  meth <- data.frame(t(feat_data[features$probe])) 
  colnames(meth) <- as.character(feat_data$SentrixID)
  meth$gene <- features$gene[match(rownames(meth),features$probe)]
  meth <- meth %>% group_by(gene)  %>% summarise(across(everything(), list(mean)))  %>% as.data.frame()
  meth <- meth[complete.cases(meth), ]
  cat(paste0("[ Genes ] : ", dim(meth)[1],"\n"))
  rownames(meth) <- meth$gene
  feat_data <- cbind(data[clincols],t(meth[2:(length(meth))]))
  return(feat_data)
}

##########
# Function to predict cancer before a given age of onset cutoff given a xgboost model
#'
#' @rdname pred_cancer_xgboost_test
#' @name pred_cancer_xgboost_test
#' 
#' @param test_dat Test data, must contain all features and outcome variable
#' @param features Methylation features in model
#' @return Dataframe containing predictions and clinical variables
##########

pred_cancer_xgboost_test <- function(test_dat, features) {

  ## Read in predictive model ##
  model <- readRDS('NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_canceratdraw_systreat_xgboost_ageofonset_model.rds')

  age_cutoff = 72
  # get intersection of features and real data
  model_features <- c('gender', 'cancer_atdraw','systemic_treatment_atdraw', features)
  model_features <- intersect(model_features, colnames(test_dat))
  test_dat$gender <- ifelse(test_dat$gender == "M", 0, 1)
  test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)

  # get y
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), "No", "Yes"))

  # get clinical data
  test_clin <- test_dat[colnames(data)[!grepl('cg',colnames(data))]]

  # get model data
  test_dat <- test_dat[, model_features]

  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'prob')

  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)

}

##########
# Function to calibrate probability scores
#'
#' @rdname platt_scaling 
#' @name platt_scaling
#' 
#' @param test_results Dataframe with uncalibrated predicted probabilites from test data
#' @return Dataframe with calibrated predicted probabilites from test data
##########

platt_scaling <- function(test_results) {
  ## Read in platt scaling recalibration model ## 
  recal_model <- readRDS('NoobCorrected_after_covs_beta_ProjPC2Adj_lfs_5UTR_scaled_canceratdraw_systreat_xgboost_recalibration_model.rds')
  # predicting on the test dataset using Platt Scaling
  ll_df_test<-data.frame(x=test_results[["test_pred.Yes"]])
  result_test_platt<-predict(recal_model,ll_df_test,type="response")
  test_results$test_pred_calibrated.Yes <- result_test_platt
  return(test_results)
}

##########
#' Function to calculate F1 score
#'
#' @rdname get_f1
#' @name get_f1
#' 
#' @param data Dataframe containing predicted probabilities and actual labels
#' @return F1 score
##########

get_f1 <- function (data) {
  precision <- Precision(data$test_label, data$predicted_label, positive = 1)
  recall <- Recall(data$test_label, data$predicted_label, positive = 1)
  f1 <- 2*precision*recall/(precision+recall)
  return(f1)
}


##########
# Function to get ROC, F1, sensitivity and specificity at optimal cutoff in test data
#'
#' @rdname ROCInfo_atcutoff
#' @name ROCInfo_atcutoff
#' 
#' @param data Dataframe containing predicted and actual values
#' @param other_title Title for plot
#' @return List containing prediction metrics and roc plot
##########

ROCInfo_atcutoff <- function(data,other_title) {

  # decision boundary 
  cutoff <- 0.5

  if (any(unique(data$test_label) %in% c("Yes","No"))) {
    data$test_label <- factor(ifelse(data$test_label == "Yes", 1, 0))
  }

  # calculate the roc values
  pred <- prediction( data$test_pred_calibrated.Yes, data$test_label )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  auc <- performance( pred, "auc" )@y.values[[1]]

  ## get sens/spec/f1
  data$predicted_label <- factor(ifelse(data$test_pred_calibrated.Yes >= cutoff,1,0))
  sensitivity <- sensitivity(data$predicted_label,data$test_label,positive=1)
  specificity <- specificity(data$predicted_label,data$test_label,positive=1)
  f1 <- get_f1(data)
   
  options(scipen = '999')

  # the main title for the plot
  sub_title <- sprintf(other_title,  "Cutoff at %.2f , AUC = %.3f", 
                       cutoff, auc )
  # plot roc curve
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) + 
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) + 
    labs( title = sub_title, x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = sensitivity, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = 1-specificity, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))

  return( list( data = data,
                plot    = roc_plot, 
                pred      = pred,
                perf      = perf,
                cutoff    = cutoff, 
                auc         = auc,
                sensitivity = sensitivity, 
                specificity = specificity,
                f1 = f1) )
}

##########
# Function that plots the top 2 PCs coloured by confounders
#'
#' @rdname generate_pcplots 
#' @name generate_pcplots 
#' 
#' @param pc_clin Dataframe with principal components (i.e. PC1, PC2, PC3, etc.) and confounder variables: array, batch, agesamplecollection 
##########

generate_pcplots <- function(pc_clin) {
  cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$gender)))
  pdf(paste0(outdir,"Plots/",output,"_PCA_gender.pdf"),width=9,height=7)
  cancerstatusplot <- ggplot(pc_clin,aes(PC1,PC2,color=cancerstatus)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(cancerstatusplot)
  suppressMessages(dev.off())

  pdf(paste0(outdir,"Plots/",output,"_PCA_agesamplecollection.pdf"),width=9,height=7)
  ageplot <- ggplot(pc_clin,aes(PC1,PC2,color=agesamplecollection)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(ageplot)
  suppressMessages(dev.off())

  cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$batch)))
  pc_clin$array <- factor(pc_clin$array)
  pdf(paste0(outdir,"Plots/",output,"_PCA_confounders.pdf"),width=9,height=7)
  confounderplot <- ggplot(pc_clin,aes(PC1,PC2,color=batch,shape=array)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(confounderplot)
  suppressMessages(dev.off())

}

