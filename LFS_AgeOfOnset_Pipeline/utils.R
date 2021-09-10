#####################################################################
############################# Libraries #############################
#####################################################################

require(dplyr)
require(reshape2)
require(minfi)
require(bumphunter)
require(dplyr)
require(minfi)
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
require(e1071)

#####################################################################
############################# Functions #############################
#####################################################################

##########
# Function to remove batch effect
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
# Function that preprocesses removes outliers
##########

remove_outliers <- function(pc,n) {
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
# Function that quantifies association between PCs and confounders/cancer status
##########

generate_pcsummary <- function(pc_clin, filename, outdir) {
  pcs <- colnames(pc_clin)[grepl( "PC" , names(pc_clin) )] ;
  anova_array <- list() ; anova_batch <- list() ; anova_cancerstatus <- list() ; anova_cancertype <- list()
  for (i in pcs) {
    anova_array[[i]] <- wilcox.test(pc_clin[,i] ~ as.factor(pc_clin$array))$p.value
    anova_batch[[i]] <- summary(aov(pc_clin[,i] ~ as.factor(pc_clin$Project)))[[1]][1,"Pr(>F)"]
    anova_cancerstatus[[i]] <- wilcox.test(pc_clin[,i] ~ as.factor(pc_clin$cancerstatus))$p.value
  }

  anova <- data.frame(PC=pcs,
                      batch_p=as.numeric(do.call("cbind",anova_batch)),
                      array_p=as.numeric(do.call("cbind",anova_array)),
                      cancerstatus_p=as.numeric(do.call("cbind",anova_cancerstatus)),
		      cancertype_p=as.numeric(do.call("cbind",anova_cancertype)))
  anova <- cbind(anova,data.frame(t(summary(pc)$importance)))
  anova$batch_padj <- anova$batch_p/anova$Proportion.of.Variance
  anova$array_padj <- anova$array_p/anova$Proportion.of.Variance
  anova$cancerstatus_padj <- anova$cancerstatus_p/anova$Proportion.of.Variance
  write.csv(anova,paste0(outdir,'Output/',filename),quote=F,row.names=F)
}

##########
# Plot top 2 PCs
##########

generate_pcplots <- function(pc_clin,output,outdir) {
  cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$cancerstatus)))
  pdf(paste0(outdir,"Plots/",output,"_PCA_cancerstatus.pdf"),width=9,height=7)
  cancerstatusplot <- ggplot(pc_clin,aes(PC1,PC2,color=cancerstatus)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(cancerstatusplot)
  suppressMessages(dev.off())

  cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$gender)))
  pdf(paste0(outdir,"Plots/",output,"_PCA_gender.pdf"),width=9,height=7)
  genderplot <- ggplot(pc_clin,aes(PC1,PC2,color=gender)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(genderplot)
  suppressMessages(dev.off())

  pdf(paste0(outdir,"Plots/",output,"_PCA_age.pdf"),width=9,height=7)
  ageplot <- ggplot(pc_clin,aes(PC1,PC2,color=agesamplecollection)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(ageplot)
  suppressMessages(dev.off())

  cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$Project)))
  pc_clin$array <- factor(pc_clin$array)
  pdf(paste0(outdir,"Plots/",output,"_PCA_confounders.pdf"),width=9,height=7)
  confounderplot <- ggplot(pc_clin,aes(PC1,PC2,color=Project,shape=array)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal()
  print(confounderplot)
  suppressMessages(dev.off())

}

##########
# Function to project 450k data onto the 850k space
##########

remove_array_confounder <- function(data) {
  cat("[ Project 450k data onto 850k space ]","\n")
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
# Function to scale methylation variables
# data = preprocessed methylation data
########## 

scale_df <- function(data,probes) {
  tmp <- data[probes]
  tmp <- scale(tmp, center=TRUE,scale=TRUE)
  data_scaled <- cbind(data[,-probes],tmp)
  return(data_scaled)
}

##########
# Function to extract probes based on location and aggregate by gene
# data = entire methylation dataset
# location = 3UTR, 5UTR, Body, TSS200, TSS1500, 1stExon
##########

aggregate_probes <- function(data,probes) {
  clincols <- colnames(data)[!grepl('cg',colnames(data))]
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- ann450k[ann450k$Name %in% probes,]
  data <- data[c(clincols,probes)]
  meth <- data.frame(t(data[probes])) 
  colnames(meth) <- as.character(data$SentrixID)
  
  meth$gene <- do.call("rbind",str_split(ann450k$UCSC_RefGene_Name[match(rownames(meth),ann450k$Name)],';'))[,1]
  meth <- meth %>% group_by(gene)  %>% summarise(across(everything(), list(mean)))  %>% as.data.frame()
  meth <- meth[complete.cases(meth), ]
  cat(paste0("[ Genes ] : ", dim(meth)[1],"\n"))
  rownames(meth) <- meth$gene
  data <- cbind(data[clincols],t(meth[2:(length(meth))]))
  return(data)
}

##########
# Function to set up parameters for bumphunter. Preprocess data to include samples of interest prior to calling this function
# data = processed methylation data
# ratio_set = ratio set objects for probes
##########

process_params <- function(data,ratio_set) {

  probes <- colnames(data)[grepl('cg',colnames(data))]
  clinical <- data[,-probes]
  #clinical <- data %>% select(ids:family_name)
  #start_loc = match("family_name",names(clinical)) + 1
  # get granges object
  g_ranges <- as.data.frame(getLocations(ratio_set))
  # get probes from rownames
  g_ranges$probe <- rownames(g_ranges)
  # remove ch and duplicatee
  g_ranges <- g_ranges[!duplicated(g_ranges$start),]
  g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]
  # beta = A matrix with rows representing genomic locations and columns representing
  beta <- t(data.frame(data[probes], row.names=data$ids))

  # beta <- t(data.frame(data[start_loc:length(data)], row.names=data$ids))
  g_ranges <- g_ranges %>%
    filter(probe %in% rownames(beta))
  beta <- beta[ rownames(beta) %in% g_ranges$probe, ]

  # pos = the genomic position
  pos <- g_ranges$start
  # chr = chromosome
  chr <- g_ranges$seqnames

  params <- list(g_ranges, beta, chr, pos, clinical)
  return(params)
}

##########
# Function to run bumphunter to identify DMRs
# params = bumphunter parameters outputted by process_params 
# designMatrix = matrix with rows representing samples and columns representing covariates. Include a column for the intercept (= 1) and a column for the variable of interest (cases vs contorls).
# name = file name
# outdir = output directory to save results to
##########

bump_hunter <- function(params, designMatrix, name, outdir) {

  DELTA_BETA_THRESH = .05 # DNAm difference threshold
  NUM_BOOTSTRAPS = 100   # number of randomizations

  # check dimensions
  stopifnot(dim(params[[2]])[2] == dim(designMatrix)[1])
  stopifnot(dim(params[[2]])[1] == length(params[[3]]))
  stopifnot(dim(params[[2]])[1] == length(params[[4]]))

  tab <- bumphunter(as.matrix(params[[2]]),
                    designMatrix,
                    chr = params[[3]],
                    pos = params[[4]],
                    nullMethod = "bootstrap",
                    cutoff = DELTA_BETA_THRESH,
                    B = NUM_BOOTSTRAPS)

  #It will provide boundaries for the differentially methylated regions, but it doesn't directly provide beta
  #values in the result (although you can get beta values for each site using other minfi functions)

  #Start and end columns indicate the limiting genomic locations of the DMR;
  #Value column indicates the average difference in methylation in the bump,
  #Area column indicates the area of the bump with respect to the 0 line.
  #Fwer column returns the family-wise error rate (FWER) of the regions estimated by the permutation scheme.
  #One can filter the results by picking a cutoff on the FWER.

  save(tab, file = paste0(outdir,'rds/',name,".rds"), compress = "xz")

}

##########
# Function to remove technical replicates 
# data = methylation data
##########

remove_duplicates <- function(data) {
  clincols <- colnames(data)[!grepl('cg',colnames(data))]
  clin <- data[clincols] %>%
    arrange(desc(array)) %>%
    distinct(tm_donor, .keep_all = TRUE) %>%
    distinct(ids, .keep_all = TRUE)
  data_cleaned <- data[data$SentrixID %in% clin$SentrixID,]
  return(data_cleaned)
}

##########
# Get ROC, F1, sensitivity and specificity at optimal cutoff in the validation set 
# data = dataframe containing predicted probabilities and actual labels
# predict = column name of predicted probabilities
# actual = column name of actual label
##########
get_auc <- function( data, predict, actual) {
  pred <- prediction( data[[predict]], data[[actual]] )
  auc <- performance( pred, "auc" )@y.values[[1]]
  return(auc)
}

##########
# Get ROC, F1, sensitivity and specificity at optimal cutoff in the validation set 
# data = dataframe containing predicted probabilities and actual labels
# predict = column name of predicted probabilities
# actual = column name of actual label
# cutoff = predicted probability cutoff
##########

get_f1 <- function (data , predict, actual, cutoff) {
  y_pred <- ifelse(data[[predict]] >= cutoff, 1, 0)
  precision <- Precision(data[[actual]], y_pred, positive = 1)
  recall <- Recall(data[[actual]], y_pred, positive = 1)
  f1 <- 2*precision*recall/(precision+recall)
  return(f1)
}


##########
# Get ROC, F1, sensitivity and specificity at optimal cutoff in the external test set 
# data = 
# predict = column name of predicted probabilities
# actual = column name of actual label
# other_title = title for plot
# model = model used i.e. xgboost, enet, rf, gbm, svm
##########

ROCInfo_atcutoff <- function( data, predict, actual, cutoff, other_title, model)
{

  if (model %in% c("rf","gbm","svmLinear","svmRadial","xgboost") ) {
    predict <- paste0(predict,".Yes")
  } 

  if (any(unique(data[[actual]]) %in% c("Yes","No"))) {
    data[[actual]] <- ifelse(data[[actual]] == "Yes", 1, 0)
  }

  # calculate the values using the ROCR require
  # true positive, false postive 
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
 
  ## get sens/spec/f1
  cm <- table(data[[actual]], ifelse(data[[predict]] >= cutoff,1,0) )
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
                plot    = plot, 
                pred      = pred,
                perf      = perf,
                cutoff    = cutoff, 
                auc         = auc,
                sensitivity = sensitivity, 
                specificity = specificity,
                f1 = f1) )
}

##########
# Function to predict cancer before a given age of onset cutoff given a xgboost model
# model = trained xgboost model 
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
##########

pred_cancer_xgboost_test <- function(model,
                                         test_dat,
                                         age_cutoff,
                                         gender,
                                         tech,
                                         cancer_atdraw,
                                         syst_treat,
                                         features) {

  # get intersection of features and real data
  intersected_feats <- intersect(features, colnames(test_dat))

  if (gender) {
    intersected_feats <- c('gender', intersected_feats)
    test_dat$gender <- ifelse(test_dat$gender == "M", 0, 1)
  }
  
  if (tech) {
    intersected_feats <- c('array', intersected_feats)
    test_dat$array <- ifelse(test_dat$array == "450", 0, 1)
  }

  if (cancer_atdraw) {
    intersected_feats <- c('cancer_atdraw', intersected_feats)
    test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  }

  if (syst_treat) {
    intersected_feats <- c('systemic_treatment_atdraw', intersected_feats)
    test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
  }

  # # get y
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), "No", "Yes"))

  # get clinical data
  test_clin <- test_dat[1:44]

  # get model data
  test_dat <- test_dat[, intersected_feats]

  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'prob')

  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)

}

##########
# Function to calibrate probability scores
# train_results = dataframe with uncalibrated predicted probabilites from training data
# test_results = dataframe with uncalibrated predicted probabilites from test data
##########

platt_scaling <- function(model_log, test_results) {
  # predicting on the test dataset using Platt Scaling
  ll_df_test<-data.frame(x=test_results[["test_pred.Yes"]])
  result_test_platt<-predict(model_log,ll_df_test,type="response")
  test_results$test_pred_calibrated.Yes <- result_test_platt
  return(test_results)
}


