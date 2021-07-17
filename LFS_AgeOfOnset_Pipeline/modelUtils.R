#####################################################################
############################# Libraries #############################
#####################################################################

require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
require(stringr)
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
# Function to scale methylation variables
# data = preprocessed methylation data
##########

scale_df <- function(data) {
  tmp <- data[45:length(data)]
  tmp <- scale(tmp, center=TRUE,scale=TRUE)
  data_scaled <- cbind(data[1:44],tmp)
  return(data_scaled)
}

##########
# Function to extract probes based on location and aggregate by gene
# data = entire methylation dataset
# location = 3UTR, 5UTR, Body, TSS200, TSS1500, 1stExon
##########

get_func_probes <- function(data,location) {
  if (location == "3UTR") {
  	location = "3'UTR"
  } 

  if (location == "5UTR") {
        location = "5'UTR"
  }
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- ann450k[ann450k$Name %in% names(data),]
  locations <- unique(do.call("rbind",str_split(ann450k$UCSC_RefGene_Group,";"))[,1])[-c(2)]
  other <- paste0(locations[!locations %in% location],collapse="|")
  ann450k <- ann450k[grepl(location,ann450k$UCSC_RefGene_Group) & !grepl(other, ann450k$UCSC_RefGene_Group),]
  cat(paste0("[ ",location, " Probes ] : ", dim(ann450k)[1],"\n")) 
  data <- data[c(colnames(data)[1:44],ann450k$Name)]
  meth <- data.frame(t(data[45:length(data)])) 
  colnames(meth) <- as.character(data$SentrixID)
  
  meth$gene <- do.call("rbind",str_split(ann450k$UCSC_RefGene_Name[match(rownames(meth),ann450k$Name)],';'))[,1]
  meth <- meth %>% group_by(gene)  %>% summarise(across(everything(), list(mean)))  %>% as.data.frame()
  meth <- meth[complete.cases(meth), ]
  cat(paste0("[ Genes ] : ", dim(meth)[1],"\n"))
  rownames(meth) <- meth$gene
  data <- cbind(data[1:44],t(meth[2:(length(meth))]))
  return(data)
}

##########
# Function to extract probes based on specified location in gene 
# data = preprocessed methylation data
# location = 3UTR, 5UTR, Body, TSS200, TSS1500, 1stExon
##########

get_func_probesonly <- function(data,location) {
  if (location == "3UTR") {
        location = "3'UTR"
  }

  if (location == "5UTR") {
        location = "5'UTR"
  }

  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- ann450k[ann450k$Name %in% names(data),]
  locations <- unique(do.call("rbind",str_split(ann450k$UCSC_RefGene_Group,";"))[,1])[-c(2)]
  other <- paste0(locations[!locations %in% location],collapse="|")
  ann450k <- ann450k[grepl(location,ann450k$UCSC_RefGene_Group) & !grepl(other, ann450k$UCSC_RefGene_Group),]
  cat(paste0("[ ",location, " Probes ] : ", dim(ann450k)[1],"\n"))
  data <- data[c(colnames(data)[1:44],ann450k$Name)]
  return(data)
}

##########
# Function to sample sample_p probes from a specified location of sample_g genes 
# data = preprocessed methylation data
# location = 3UTR, 5UTR, Body, TSS200, TSS1500, 1stExon
# sample_g = number of genes
# sample_p = number of probes
# seed = seed # 
##########

sample_func_probes <- function(data,location,sample_g,sample_p,seed) {
  if (location == "3UTR") {
        location = "3'UTR"
  }

  if (location == "5UTR") {
        location = "5'UTR"
  }
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ## extract all locations 
  ann450k <- ann450k[ann450k$Name %in% names(data),]
# locations <- c("Body","TSS200","1stExon","TSS1500","5'UTR","3'UTR")
  locations <- unique(do.call("rbind",str_split(ann450k$UCSC_RefGene_Group,";"))[,1])[-c(2)]
  other <- paste0(locations[!locations %in% location],collapse="|")
  ## get all locations excluding the one of interest
  ann450k <- ann450k[grepl(location,ann450k$UCSC_RefGene_Group) & !grepl(other, ann450k$UCSC_RefGene_Group),]
  cat(paste0("[ ",location," Probes ] :",dim(ann450k)[1],"\n"))
  set.seed(seed)
  ## get all probes in functional region
  data <- data[c(colnames(data)[1:44],ann450k$Name)]
  meth <- data.frame(t(data[45:length(data)]))
  colnames(meth) <- as.character(data$SentrixID)
  meth$gene <- do.call("rbind",str_split(ann450k$UCSC_RefGene_Name[match(rownames(meth),ann450k$Name)],';'))[,1]
  ## sample number of genes
  genes <-  sample(unique(meth$gene),sample_g,replace=FALSE)
  cat(paste0("[ Sampled Genes ] :",length(genes),"\n"))
  meth <- meth[meth$gene %in% genes,]
  meth$probe <- rownames(meth)
  ## sample number of probes
  sampled <- meth %>% group_by(gene) %>% sample_n(1) %>% data.frame()
  samp_probes1 <- sampled$probe
  sampled <- meth[!meth$probe %in% samp_probes1,] %>% group_by(gene) %>% sample_n(1) %>% data.frame()
  remain_p <- sample_p - sample_g
  if (length(sampled$probe) - remain_p > 0 ) {
     samp_probes2  <- sample(sampled$probe,remain_p,replace=FALSE)
  } else {
     samp_probes2  <- sampled$probe
  }
  samp_probes <- c(samp_probes1,samp_probes2)
  cat(paste0("[ Sampled Probes ] :",length(samp_probes),"\n"))
  sampled_data <- data[c(colnames(data)[1:44],samp_probes)]
  return(sampled_data)
}

##########
# Function to aggregate sampled probes by gene for sampling analysis
# data = preprocessed methylation data
##########

aggregate_sampled_probes <- function(data) {
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- ann450k[ann450k$Name %in% names(data),]
  cat(paste0("[ Probes ] : ", dim(ann450k)[1],"\n"))
  data <- data[c(colnames(data)[1:44],ann450k$Name)]
  meth <- data.frame(t(data[45:length(data)]))
  colnames(meth) <- as.character(data$SentrixID)
  meth$gene <- do.call("rbind",str_split(ann450k$UCSC_RefGene_Name[match(rownames(meth),ann450k$Name)],';'))[,1]
  meth <- meth %>% group_by(gene)  %>% summarise(across(everything(), list(mean)))  %>% as.data.frame()
  cat(paste0("[ Genes ] : ", dim(meth)[1],"\n"))
  rownames(meth) <- meth$gene
  data <- cbind(data[1:44],t(meth[2:(length(meth))]))
  return(data)
}

##########
# Function to get probes differentially methylated between wildtype and mutant patients
# data = preprocessed methylation data
##########

get_lfs_probes <- function(data) {
  lfsprobes <- readRDS("Resources/lfs_probes.rds")
  keep <- c(colnames(data)[1:44],lfsprobes[lfsprobes %in% colnames(data)])
  data_lfs <- data[keep]
  return(data_lfs)
}

##########
# Function to get probes within gene list formed by Daniel Cole
# data = preprocessed methylation data
##########

get_gene_probes <- function(data) {
  geneprobes <- readRDS("Resources/gene_probes.rds")
  keep <- c(colnames(data)[1:44],geneprobes[geneprobes %in% colnames(data)])
  data_gene <- data[keep]
  return(data_gene)
}

##########
# Function to remove cancer signal using bumphunter
# data = preprocessed methylation data
# id = name of file
# outdir = output directory to save results
##########

get_cancer_probes <- function(data, id, outdir) {
  registerDoParallel(cores = 4)

  # the desgin matrix =  with rows representing samples and columns representing covariates.
  # regression is applied to each row of mat. basically a column for the intercept (= 1)
  # and a column for the variable of interest (cases vs controls).
  designMatrix <- data.frame(intcpt = 1, cancer = ifelse(data$cancerstatus == "Unaffected",0,1))

  # get granges object
  ratio_set <- readRDS('Resources/raw_ratio_set.rds')
  cat(paste0("[ Process parameters ]","\n"))
  params <- process_params(data,ratio_set)

  if(!file.exists(paste0(outdir,'rds/',id,"bumps.rds"))){
	cat(paste0("[ Generating cancer probes ]","\n"))  
	bumps <- bump_hunter(params, designMatrix, paste0(id,"bumps"),outdir)
  } else {
        cat(paste0("[ Loading cancer probes ]","\n"))
  	bumps <- get(load(paste0(outdir,'rds/',id,"bumps.rds")))
  }
  g_ranges <- as.data.frame(getLocations(ratio_set))
  # get probes from rownames
  g_ranges$probe <- rownames(g_ranges)
  # remove ch and duplicatee
  g_ranges <- g_ranges[!duplicated(g_ranges$start),]

  #bumps <- get(load("Data_objects/bumps.rda"))
  bumps_table <- bumps$table
  bumps_sig <- bumps_table[bumps_table$p.value < 0.05,]
  cat(paste0("Number of bumps : ",dim(bumps_sig)[1],"\n"))
  bumps_probes <- inner_join(bumps_sig, g_ranges)$probe
  return(bumps_probes)
}

##########
# Function to set up parameters for bumphunter. Preprocess data to include samples of interest prior to calling this function
# data = processed methylation data
# ratio_set = ratio set objects for probes
##########

process_params <- function(data,ratio_set) {

  clinical <- data[1:44]
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
  beta <- t(data.frame(data[45:length(data)], row.names=data$ids))

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
# Function to remove cross-reactive probes
# data = methylation data
##########

remove_xRP <- function(data) {
  xReactiveProbes <- read.csv(file="Resources/cross_reactive_probes.csv", stringsAsFactors=FALSE)
  data_NoxRP <- data[colnames(data)[!colnames(data) %in% xReactiveProbes$TargetID]]
  return(data_NoxRP)
}

##########
# Function to remove age-associated probes 
# data = methylation data
##########

remove_age_probes <- function(data) {
  age_probes <- read.csv('Resources/Horvath_PedBE_probes.csv', stringsAsFactors=FALSE)
  data <- data[colnames(data)[!colnames(data) %in% age_probes$probe]]
  return(data)
}

##########
# Function to remove technical replicates 
# data = methylation data
##########

remove_duplicates <- function(data) {
  clin <- data[1:44] %>%
  	arrange(desc(array)) %>%
  	distinct(tm_donor, .keep_all = TRUE) %>%
  	distinct(ids, .keep_all = TRUE)
  data_cleaned <- data[data$SentrixID %in% clin$SentrixID,]
  return(data_cleaned)
}

##########
# Function to build an elastic net model to predict cancer before a given age of onset cutoff  
# train_dat = training data
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
##########

pred_cancer_enet <- function(train_dat, 
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			                       cancer_atdraw,
			                       syst_treat,
                             features) {
  
  
  # get intersection of bh features and real data
  features <- as.character(unlist(features))
  
  intersected_feats <- intersect(features, colnames(train_dat))
  
  if(gender) {
    intersected_feats <- c('gender', intersected_feats)
    train_dat$gender <- ifelse(train_dat$gender == "M", 0, 1)
    test_dat$gender <- ifelse(test_dat$gender == "M", 0, 1)
  }
  if (tech) {
    intersected_feats <- c('array', intersected_feats)
    train_dat$array <- ifelse(train_dat$array == "450", 0, 1)
    test_dat$array <- ifelse(test_dat$array == "450", 0, 1)
  }
  if (cancer_atdraw) {
    intersected_feats <- c('cancer_atdraw', intersected_feats)
    train_dat$cancer_atdraw <- ifelse(train_dat$cancer_atdraw == 'Yes', 1, 0)
    test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  }

  if (syst_treat) {
    intersected_feats <- c('systemic_treatment_atdraw', intersected_feats)
    train_dat$systemic_treatment_atdraw <- ifelse(train_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
    test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
  }
  
  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
  # # get y
  train_y <- factor(ifelse(train_dat$ageofonset > age_cutoff | is.na(train_dat$ageofonset), 0, 1))
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), 0, 1))

  # get clinical data
  train_clin <- train_dat[1:44]
  test_clin <- test_dat[1:44]
 
  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  N_CV_REPEATS = 2
  nfolds = 5
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(train_dat)
                                                , y =  train_y
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = type_measure
                                                , family = type_family
                                                , standardize = FALSE 
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
      )
      elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
    # as your value for alpha
  }
  
  # stop if you did not recover error for any models 
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  
  # get index of best alpha (lowest error) - alpha is values 0.1-0.9
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(train_dat)
                                     , y =  train_y
                                     , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    # get optimal lambda - the tuning parameter for ridge and lasso
    # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
    # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
    # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
    # GIVE YOU REASONS
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(train_dat)
                  , y =  train_y
                  ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  
  return(list(model, temp_dat, best_alpha, temp.min_lambda_index))  
  
}

##########
# Function to build xgboost model to predict cancer before a given age of onset cutoff 
# train_dat = training data
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
##########

pred_cancer_xgboost <- function(train_dat,
                                test_dat,
                                age_cutoff,
                                gender,
                                tech,
                                cancer_atdraw,
                                syst_treat,
                                features) {


  # get intersection of bh features and real data
  intersected_feats <- intersect(features, colnames(train_dat))

  if(gender) {
    intersected_feats <- c('gender', intersected_feats)
    train_dat$gender <- ifelse(train_dat$gender == "M", 0, 1)
    test_dat$gender <- ifelse(test_dat$gender == "M", 0, 1)
  }
  if (tech) {
    intersected_feats <- c('array', intersected_feats)
    train_dat$array <- ifelse(train_dat$array == "450", 0, 1)
    test_dat$array <- ifelse(test_dat$array == "450", 0, 1)
  }

  if (cancer_atdraw) {
    intersected_feats <- c('cancer_atdraw', intersected_feats)
    train_dat$cancer_atdraw <- ifelse(train_dat$cancer_atdraw == 'Yes', 1, 0)
    test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  }

  if (syst_treat) {
    intersected_feats <- c('systemic_treatment_atdraw', intersected_feats)
    train_dat$systemic_treatment_atdraw <- ifelse(train_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
    test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
  }

  train_y <- factor(ifelse(train_dat$ageofonset > age_cutoff | is.na(train_dat$ageofonset), "No", "Yes"))
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), "No", "Yes"))

  # get clinical data
  train_clin <- train_dat[1:44]
  test_clin <- test_dat[1:44]

  train_dat <- train_dat[intersected_feats]
  test_dat <- test_dat[intersected_feats]

  # determines how you train the model.
  NFOLDS <- 5
  fitControl <- trainControl(
    method = "repeatedcv",
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary

  )

  nrounds = 1000

  gbmGrid <- expand.grid(
    nrounds = seq(from = 200, to = nrounds, by = 200),
    eta = c(0.025, 0.05, 0.1, 0.3),
    max_depth = c(2, 4, 6),
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = 1,
    subsample = 1
  )

  gbm.model <- train(x = train_dat
                     , y = train_y
                     , method = "xgbTree"
                     , metric = 'ROC'
                     , tuneGrid = gbmGrid
                     , trControl = fitControl
                     , verbose = FALSE
  )

  test.predictions <- predict(gbm.model
                              ,newdata = test_dat
                              ,type='prob')

  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))

  return(list(gbm.model,temp_dat))

}

##########
# Function to build model to predict cancer before a given age of onset cutoff using a gbm
# train_dat = training data
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
##########

pred_cancer_gbm <- function(train_dat,
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
                             cancer_atdraw,
                             syst_treat,
                             features) {


  # get intersection of bh features and real data
  features <- as.character(unlist(features))

  intersected_feats <- intersect(features, colnames(train_dat))

  if(gender) {
    intersected_feats <- c('gender', intersected_feats)
    train_dat$gender <- ifelse(train_dat$gender == "M", 0, 1)
    test_dat$gender <- ifelse(test_dat$gender == "M", 0, 1)
  }
  if (tech) {
    intersected_feats <- c('array', intersected_feats)
    train_dat$array <- ifelse(train_dat$array == "450", 0, 1)
    test_dat$array <- ifelse(test_dat$array == "450", 0, 1)
  }
  if (cancer_atdraw) {
    intersected_feats <- c('cancer_atdraw', intersected_feats)
    train_dat$cancer_atdraw <- ifelse(train_dat$cancer_atdraw == 'Yes', 1, 0)
    test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  }

  if (syst_treat) {
    intersected_feats <- c('systemic_treatment_atdraw', intersected_feats)
    train_dat$systemic_treatment_atdraw <- ifelse(train_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
    test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
  }
  # # get y
  train_y <- factor(ifelse(train_dat$ageofonset > age_cutoff | is.na(train_dat$ageofonset), "No", "Yes"))
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), "No", "Yes"))

  # get clinical data
  train_clin <- train_dat[1:44]
  test_clin <- test_dat[1:44]

  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]

  # determines how you train the model.
  NFOLDS <- 5
  fitControl <- trainControl(
    method = "repeatedcv",
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary

  )

  gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 6, 9),
                    n.trees = (0:25)*50, 
                    shrinkage = c(0.001, 0.01, 0.1),
                    n.minobsinnode = 10) 

  gbm.model <- train(x = train_dat
                           , y = train_y
                           , method = "gbm"
                           , metric = 'ROC'
			   , tuneGrid = gbmGrid
                           , trControl = fitControl
                           , verbose = FALSE
  )

  test.predictions <- predict(gbm.model
                              ,newdata = test_dat
                              ,type='prob')

  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))

  return(list(gbm.model,temp_dat))

}

##########
# Function to build svm model to predict cancer before a given age of onset cutoff
# train_dat = training data
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
##########

pred_cancer_svm <- function(train_dat,
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			                       cancer_atdraw,
			                       syst_treat,
                             features,
			                       kernel) {

  # get intersection of bh features and real data
  features <- as.character(unlist(features))

  intersected_feats <- intersect(features, colnames(train_dat))

  if(gender) {
    intersected_feats <- c('gender', intersected_feats)
    train_dat$gender <- ifelse(train_dat$gender == "M", 0, 1)
    test_dat$gender <- ifelse(test_dat$gender == "M", 0, 1)
  }
  if (tech) {
    intersected_feats <- c('array', intersected_feats)
    train_dat$array <- ifelse(train_dat$array == "450", 0, 1)
    test_dat$array <- ifelse(test_dat$array == "450", 0, 1)
  }
  if (cancer_atdraw) {
    intersected_feats <- c('cancer_atdraw', intersected_feats)
    train_dat$cancer_atdraw <- ifelse(train_dat$cancer_atdraw == 'Yes', 1, 0)
    test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  }

  if (syst_treat) {
    intersected_feats <- c('systemic_treatment_atdraw', intersected_feats)
    train_dat$systemic_treatment_atdraw <- ifelse(train_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
    test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
  }

  # get y
  train_y <- factor(ifelse(train_dat$ageofonset > age_cutoff | is.na(train_dat$ageofonset), 0, 1))
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), 0, 1))

  # get clinical data
  train_clin <- train_dat[1:44]
  test_clin <- test_dat[1:44]

  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]

  # determines how you train the model.
  NFOLDS <- 5
  fitControl <- trainControl(
    method = "repeatedcv",  
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary

  )

  svm.model <- train(x = train_dat
                           , y = train_y
                           , method = kernel
			                     , metric = 'ROC'
                           , trControl = fitControl
                           , verbose = FALSE                       
  )
  
  test.predictions <- predict(svm.model 
			      ,newdata = test_dat
			      ,type='prob')

  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))

  return(list(svm.model,temp_dat))

}

##########
# Function to build rf model to predict cancer before a given age of onset cutoff
# train_dat = train_full
# test_dat = test_full
# age_cutoff = 72
# gender = gender
# tech = tech
# features = remaining_features
##########

pred_cancer_rf <- function(train_dat, 
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			                       cancer_atdraw,
			                       syst_treat,
                             features) {
  
  
  # get intersection of bh features and real data
  features <- as.character(unlist(features))
  
  intersected_feats <- intersect(features, colnames(train_dat))
  
  if(gender) {
    intersected_feats <- c('gender', intersected_feats)
    train_dat$gender <- factor(ifelse(train_dat$gender == "M", 0, 1))
    test_dat$gender <- factor(ifelse(test_dat$gender == "M", 0, 1))
  }
  if (tech) {
    intersected_feats <- c('array', intersected_feats)
    train_dat$array <- factor(ifelse(train_dat$array == "450", 0, 1))
    test_dat$array <- factor(ifelse(test_dat$array == "450", 0, 1))
  }
  if (cancer_atdraw) {
    intersected_feats <- c('cancer_atdraw', intersected_feats)
    train_dat$cancer_atdraw <- ifelse(train_dat$cancer_atdraw == 'Yes', 1, 0)
    test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  }

  if (syst_treat) {
    intersected_feats <- c('systemic_treatment_atdraw', intersected_feats)
    train_dat$systemic_treatment_atdraw <- ifelse(train_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
    test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
  }
  
  # # get y
  train_y <- factor(ifelse(train_dat$ageofonset > age_cutoff | is.na(train_dat$ageofonset), "No", "Yes"))
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), "No", "Yes"))    

  # get clinical data
  train_clin <- train_dat[1:44]
  test_clin <- test_dat[1:44]
  
  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]

  # determines how you train the model.
  NFOLDS <- 5
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary
    
  )
  
  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  
  mtry <- sqrt(ncol(train_dat[,colnames(train_dat)]))
  tunegrid <- expand.grid(.mtry=mtry)#, .ntree=c(500,1000,1500,2000),.min.node.size = c(1, 3, 5))
  
  model <- train(x = train_dat
                 , y = train_y
                 , metric = 'ROC'
                 , method = "rf"
                 , trControl = fitControl
                 , tuneGrid = tunegrid
                 , importance = T
                 , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$X1)
  
  # Predictions on test data
  
  # This returns 100 prediction with 1-100 lambdas
  test.predictions <- predict(model, 
                              data.matrix(test_dat),
                              type = 'prob')
  
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  
  return(list(model, temp_dat, importance))
   
}

##########
# Function to predict cancer before a given age of onset cutoff given a gbm model
# model = trained gbm model 
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
##########

pred_cancer_gbm_test <- function(model,
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			                       cancer_atdraw,
			                       syst_treat,
                             features) {

  # get intersection of bh features and real data
  features <- as.character(unlist(features))

  intersected_feats <- intersect(features, colnames(test_dat))

  if(gender) {
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
# Function to predict cancer before a given age of onset cutoff given a svm model
# model = trained svm model 
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
# kernel = linear, radial, polynomial, sigmoid
##########

pred_cancer_svm_test <- function(model,
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			                       cancer_atdraw,
			                       syst_treat,
                             features, 
			                       kernel) {

  # get intersection of bh features and real data
  features <- as.character(unlist(features))

  intersected_feats <- intersect(features, colnames(test_dat))

  if(gender) {
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
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), 0, 1))

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
# Function to predict cancer before a given age of onset cutoff given a random forest model
# model = trained random forest model 
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
##########

pred_cancer_rf_test <- function(model,
			                       test_dat,
			                       age_cutoff,
                             gender,
                             tech,
			                       cancer_atdraw,
			                       syst_treat,
                             features) {

  features <- as.character(unlist(features))
  # get intersection of bh features and real data
  intersected_feats <- intersect(features, colnames(test_dat))

  if(gender) {
    intersected_feats <- c('gender', intersected_feats)
    test_dat$gender <- factor(test_dat$gender)
  }
  if (tech) {
    intersected_feats <- c('array', intersected_feats)
    test_dat$array <- factor(test_dat$array)
  }

  if (cancer_atdraw) {
    intersected_feats <- c('cancer_atdraw', intersected_feats)
    test_dat$cancer_atdraw <- ifelse(test_dat$cancer_atdraw == 'Yes', 1, 0)
  }

  if (syst_treat) {
    intersected_feats <- c('systemic_treatment_atdraw', intersected_feats)
    test_dat$systemic_treatment_atdraw <- ifelse(test_dat$systemic_treatment_atdraw == 'Yes', 1, 0)
  }

  # intersected_feats_rand <- intersect(rand_feats, colnames(training_dat))
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
# Function to predict cancer before a given age of onset cutoff given a elastic net model
# model = trained elastic net model 
# lambda_index = optimal lambda determined using CV 
# test_dat = test data
# age_cutoff = 72
# gender = 0 for M, 1 for F
# tech = 0 for 450k, 1 for 850k
# cancer_atdraw = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# syst_treat = 1 if cancer developed prior to/at draw and 0 if cancer developed after
# features = features of interest
##########

pred_cancer_enet_test <- function(model,
			                       lambda_index,
                             test_dat,
			                       age_cutoff,
                             gender,
                             tech,
			                       cancer_atdraw,
			                       syst_treat,
                             features) {

  features <- as.character(unlist(features))

  # get intersection of bh features and real data
  intersected_feats <- intersect(features, colnames(test_dat))
  
  if(gender) {
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
  test_y <- factor(ifelse(test_dat$ageofonset > age_cutoff | is.na(test_dat$ageofonset), 0, 1))

  # get clinical data
  test_clin <- test_dat[1:44]

  # get model data
  test_dat <- test_dat[, intersected_feats]

  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'response')

  test.predictions <- test.predictions[,lambda_index]

  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)
}


opt.cut = function(pred,perf){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
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
# Get ROC, F1, sensitivity and specificity at optimal cutoff in the validation set 
# data = 
# predict = column name of predicted probabilities
# actual = column name of actual label
# cost.fp = weight to represent cost of a false positive
# cost.fn = weight to represent cost of a false negative
# other_title = title for plot
# model = model used i.e. xgboost, enet, rf, gbm, svm
##########

ROCInfo_atopt <- function(data, predict, actual, cost.fp, cost.fn, other_title, model)
{
  
  if (model %in% c("rf","gbm","svmLinear","svmRadial","xgboost") ) {
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
		plot 	  = plot, 
                pred      = pred,
                perf      = perf,
                cutoff 	  = best_cutoff, 
                auc         = auc,
                sensitivity = best_tpr, 
                specificity = 1 - best_fpr,
                f1 = f1,
                totalcost   = best_cost) )
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
		            plot 	  = plot, 
                pred      = pred,
                perf      = perf,
                cutoff 	  = cutoff, 
                auc         = auc,
                sensitivity = sensitivity, 
                specificity = specificity,
                f1 = f1) )
}

##########
# Gather all metrics across the different models and aggregate into a dataframe 
# dir = directory containing all prediction metric objects
# ext = extension indicating which files to obtain (i.e. ROCInfoVal.rds for validation or ROCInfoTest.rds for test)
##########

get_all_auc <- function(dir,ext) {
  auc_enet <- list() ; sens_enet <- list() ; spec_enet <- list() ; f1_enet <- list() 
  all_models <- list.files(dir,ext)
  for (model in all_models) {
    cat(paste0(model,'\n'))
    r_all <- readRDS(paste0(dir,model))
    id <- gsub("_ageofonset_ROCInfoTest.rds|_ageofonset_ROCInfoVal.rds","",model)
    auc_enet[id] <- r_all[[6]]
    sens_enet[id] <- r_all[[7]]
    spec_enet[id] <- r_all[[8]]
    f1_enet[id] <- r_all[[9]]

  }
  auc_all <- data.frame(auc = do.call(rbind,auc_enet),
			 sens = do.call(rbind,sens_enet),
			 spec = do.call(rbind,spec_enet),
			 f1 = do.call(rbind,f1_enet))
  auc_all$preprocessing <- rownames(auc_all)
  return(auc_all)
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

platt_scaling <- function(train_results, test_results,model) {
  actual <- "test_label"
  if (model %in% c("rf","svmLinear","xgboost","svmRadial","gbm") ) {
    predict <- "test_pred.Yes"
  } else {
    predict <- "test_pred"
  }
  # calculating Log Loss without Platt Scaling
  ll <- LogLoss(train_results[[predict]],train_results[[actual]])
  ll_df <- data.frame(x=train_results[[predict]],y=as.factor(train_results[[actual]]))
  model_log <-glm(y~x,data = ll_df,family = binomial)

  # predicting on the cross validation after platt scaling
  result_platt<-predict(model_log,ll_df["x"],type = "response")
  train_results$test_pred_calibrated.Yes <- result_platt
  
  # predicting on the test dataset using Platt Scaling
  ll_df_test<-data.frame(x=test_results[[predict]])
  result_test_platt<-predict(model_log,ll_df_test,type="response")
  test_results$test_pred_calibrated.Yes <- result_test_platt
  
  return(list(train_results,test_results))
}
