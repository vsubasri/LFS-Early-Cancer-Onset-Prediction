suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library(stringr))
suppressMessages(library(bumphunter))
suppressMessages(library(dplyr))
suppressMessages(library(minfi))
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
suppressMessages(library(doParallel))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ROCR))
suppressMessages(library(data.table))
suppressMessages(library(caret))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(MLmetrics))
suppressMessages(library(e1071))
suppressMessages(library(PRROC))

# Define F1 Score function
f1Summary <- function(data, lev = NULL, model = NULL) {
  precision <- posPredValue(data$pred, data$obs, positive = lev[1])
  recall <- sensitivity(data$pred, data$obs, positive = lev[1])
  f1 <- 2*((precision*recall) / (precision+recall))
  out <- c(F1 = f1, Precision = precision, Recall = recall)
  names(out) <- c("F1", "Precision", "Recall")
  out
}

auprcSummary <- function(data, lev = NULL, model = NULL){
  the_curve <- pr.curve(data$pred, data$obs, curve=FALSE)
  out <- the_curve$auc.integral
  names(out) <- "AUPRC"
  out
  
}

scale_df <- function(data) {
  ind <- grep("array", colnames(test_dat))[1]
  tmp <- data[(ind+1):length(data)]
  tmp <- scale(tmp, center=TRUE,scale=TRUE)
  data_scaled <- cbind(data[1:ind],tmp)
  return(data_scaled)
}

##########
# function to extract probes based on location and aggregate by gene
##########
get_func_probes <- function(data,location) {
  if (location == "3UTR") {
  	location = "3'UTR"
  } 

  if (location == "5UTR") {
        location = "5'UTR"
  }

  ind <- grep("cg", colnames(data))[1]
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- ann450k[ann450k$Name %in% names(data),]
#  locations <- c("Body","TSS200","1stExon","TSS1500","5'UTR","3'UTR")
  locations <- unique(do.call("rbind",str_split(ann450k$UCSC_RefGene_Group,";"))[,1])[-c(2)]
  other <- paste0(locations[!locations %in% location],collapse="|")
  ann450k <- ann450k[grepl(location,ann450k$UCSC_RefGene_Group) & !grepl(other, ann450k$UCSC_RefGene_Group),]
  cat(paste0("[ ",location, " Probes ] : ", dim(ann450k)[1],"\n")) 
  data <- data[c(colnames(data)[1:(ind-1)],ann450k$Name)]
  meth <- data.frame(t(data[ind:length(data)])) 
  colnames(meth) <- as.character(data$SentrixID)
  
  meth$gene <- do.call("rbind",str_split(ann450k$UCSC_RefGene_Name[match(rownames(meth),ann450k$Name)],';'))[,1]
  meth <- meth %>% group_by(gene)  %>% summarise(across(everything(), list(mean)))  %>% as.data.frame()
  meth <- meth[complete.cases(meth), ]
  cat(paste0("[ Genes ] : ", dim(meth)[1],"\n"))
  rownames(meth) <- meth$gene
  data <- cbind(data[1:(ind-1)],t(meth[2:(length(meth))]))
  return(data)
}

##########
# function to extract probes based on location
##########
get_func_probesonly <- function(data,location) {
  if (location == "3UTR") {
        location = "3'UTR"
  }

  if (location == "5UTR") {
        location = "5'UTR"
  }

  ind <- grep("cg", colnames(data))[1]
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- ann450k[ann450k$Name %in% names(data),]
#  locations <- c("Body","TSS200","1stExon","TSS1500","5'UTR","3'UTR")
  locations <- unique(do.call("rbind",str_split(ann450k$UCSC_RefGene_Group,";"))[,1])[-c(2)]
  other <- paste0(locations[!locations %in% location],collapse="|")
  ann450k <- ann450k[grepl(location,ann450k$UCSC_RefGene_Group) & !grepl(other, ann450k$UCSC_RefGene_Group),]
  cat(paste0("[ ",location, " Probes ] : ", dim(ann450k)[1],"\n"))
  data <- data[c(colnames(data)[1:(ind-1)],ann450k$Name)]
  return(data)
}

sample_func_probes <- function(data,location,sample_g,sample_p,seed) {
  if (location == "3UTR") {
        location = "3'UTR"
  }

  if (location == "5UTR") {
        location = "5'UTR"
  }

  ind <- grep("cg", colnames(data))[1]-1
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ## extract all locations 
  ann450k <- ann450k[ann450k$Name %in% names(data),]
# locations <- c("Body","TSS200","1stExon","TSS1500","5'UTR","3'UTR")
  if (location != "lfs") {
  	locations <- unique(do.call("rbind",str_split(ann450k$UCSC_RefGene_Group,";"))[,1])[-c(2)]
 	 other <- paste0(locations[!locations %in% location],collapse="|")
  ## get all locations excluding the one of interest
  	ann450k <- ann450k[grepl(location,ann450k$UCSC_RefGene_Group) & !grepl(other, ann450k$UCSC_RefGene_Group),]
  }
  cat(paste0("[ ",location," Probes ] :",dim(ann450k)[1],"\n"))
  set.seed(seed)
  ## get all probes in functional region
  data <- data[c(colnames(data)[1:ind],ann450k$Name)]
  meth <- data.frame(t(data[(ind+1):length(data)]))
  colnames(meth) <- as.character(data$SentrixID)
  meth$gene <- do.call("rbind",str_split(ann450k$UCSC_RefGene_Name[match(rownames(meth),ann450k$Name)],';'))[,1]
  ## sample number of genes
  print(paste0("Number of probes available to sample: ",length(unique(meth$gene))))
  REPLACE_FLAG = FALSE ## set false to sample without replacement
  genes <-  sample(unique(meth$gene),sample_g,replace=REPLACE_FLAG)
  cat(paste0("[ Sampled Genes ] :",length(genes),"\n"))
  meth <- meth[meth$gene %in% genes,]
  meth$probe <- rownames(meth)
  ## sample number of probes 
  cat(paste0("[ Sampled Probes ] :",length(genes),"\n"))
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
  sampled_data <- data[c(colnames(data)[1:ind],samp_probes)]
  return(sampled_data)
}

get_probe_genes <- function(data) {
  ind <- grep("cg", colnames(data))[1]-1
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- ann450k[ann450k$Name %in% names(data),]
  cat(paste0("[ Probes ] : ", dim(ann450k)[1],"\n"))
  data <- data[c(colnames(data)[1:ind],ann450k$Name)]
  meth <- data.frame(t(data[(ind+1):length(data)]))
  colnames(meth) <- as.character(data$SentrixID)
  meth$gene <- do.call("rbind",str_split(ann450k$UCSC_RefGene_Name[match(rownames(meth),ann450k$Name)],';'))[,1]
  meth <- meth[c("gene","Name")]
  colnames(meth) <- c("gene", "probe")
  write.table(meth,'features.txt',sep='\t')
  return(meth)
}

aggregate_sampled_probes <- function(data) {
  ind <- grep("cg", colnames(data))[1]-1
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- ann450k[ann450k$Name %in% names(data),]
  cat(paste0("[ Probes ] : ", dim(ann450k)[1],"\n"))
  data <- data[c(colnames(data)[1:ind],ann450k$Name)]
  meth <- data.frame(t(data[(ind+1):length(data)]))
  colnames(meth) <- as.character(data$SentrixID)
  meth$gene <- do.call("rbind",str_split(ann450k$UCSC_RefGene_Name[match(rownames(meth),ann450k$Name)],';'))[,1]
  meth <- meth %>% group_by(gene)  %>% summarise(across(everything(), list(mean)))  %>% as.data.frame()
  cat(paste0("[ Genes ] : ", dim(meth)[1],"\n"))
  rownames(meth) <- meth$gene
  data <- cbind(data[1:ind],t(meth[2:(length(meth))]))
  return(data)
}

##########
# get lfs probes
##########
get_lfs_probes <- function(data) {
  ind <- grep("cg", colnames(data))[1]-1
  lfsprobes <- readRDS("Resources/lfs_probes.rds")
  keep <- c(colnames(data)[1:ind],lfsprobes[lfsprobes %in% colnames(data)])
  data_lfs <- data[keep]
  return(data_lfs)
}

##########
# get gene probes
##########
get_gene_probes <- function(data) {
  ind <- grep("cg", colnames(data))[1]-1
  geneprobes <- readRDS("Resources/gene_probes.rds")
  keep <- c(colnames(data)[1:ind],geneprobes[geneprobes %in% colnames(data)])
  data_gene <- data[keep]
  return(data_gene)
}

##########
# remove cancer signal
##########
get_cancer_probes <- function(data, id, outdir) {
  ind <- grep("cg", colnames(data))[1]-1
  registerDoParallel(cores = 4)

  # the desgin matrix =  with rows representing samples and columns representing covariates.
  # regression is applied to each row of mat. basically a column for the intercept (= 1)
  # and a column for the variable of interest (cases vs controls).
  designMatrix <- data.frame(intcpt = 1,
                           cancer = ifelse(data$cancerstatus == "Unaffected",0,1))

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
# process parameters for bumphunter
##########
process_params <- function(data,ratio_set) {
  ind <- grep("cg", colnames(data))[1]-1
  clinical <- data[1:ind]
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
  beta <- t(data.frame(data[(ind+1):length(data)], row.names=data$ids))

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
# run bumphunter
##########
bump_hunter <- function(params, designMatrix, name, outdir) {
  DELTA_BETA_THRESH = .1 # DNAm difference threshold
  NUM_BOOTSTRAPS = 100   # number of randomizations

  # check dimensions
  stopifnot(dim(params[[2]])[2] == dim(designMatrix)[1])
  stopifnot(dim(params[[2]])[1] == length(params[[3]]))
  stopifnot(dim(params[[2]])[1] == length(params[[4]]))

  # beta = A matrix with rows representing genomic locations and columns representing
  # the desgin matrix =  with rows representing samples and columns representing covariates.
  # regression is applied to each row of mat. basically a column for the intercept (= 1)
  # and a column for th variable of interest (cases vs contorls).
  # pos = the genomic position
  # chr = chromosome

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

##########b
# remove cross-reactive probes
##########

remove_xRP <- function(data) {
  xReactiveProbes <- read.csv(file="Resources/cross_reactive_probes.csv", stringsAsFactors=FALSE)
  data_NoxRP <- data[colnames(data)[!colnames(data) %in% xReactiveProbes$TargetID]]
  return(data_NoxRP)
}

##########b
# remove age-related probes 
##########

remove_age_probes <- function(data) {
  age_probes <- read.csv('Resources/Horvath_PedBE_probes.csv', stringsAsFactors=FALSE)
  data <- data[colnames(data)[!colnames(data) %in% age_probes$probe]]
  return(data)
}

remove_duplicates <- function(data) {
  ind <- grep("cg", colnames(data))[1]-1 
  clin <- data[1:ind] %>%
  	arrange(desc(array)) %>%
  	distinct(tm_donor, .keep_all = TRUE) %>%
  	distinct(ids, .keep_all = TRUE)
  data_cleaned <- data[data$SentrixID %in% clin$SentrixID,]
  return(data_cleaned)
}

get_cell_props <- function(data) {
  cell_props <- readRDS('Resources/idol_cell_props.rds')
  cell_props <- cell_props[c("ids","NK","CD8T","CD4T","Bcell","Mono","Neu")]
  data <- merge(data,cell_props, by="ids",all.x=T)
  return(data)
}

##########
# get features
##########
get_features <- function(gender,
			 tech,
			 cancer_atdraw,
			 syst_treat,
			 cell_props,
			 family,
			 tp53,
			 features) {

  intersected_feats <- as.character(unlist(features))

  if(gender) {
    intersected_feats <- c('gender', intersected_feats)

  }
  if (tech) {
    intersected_feats <- c('array', intersected_feats)
  }
  if (cancer_atdraw) {
    intersected_feats <- c('cancer_atdraw', intersected_feats)
  }

  if (syst_treat) {
    intersected_feats <- c('systemic_treatment_atdraw', intersected_feats)
  }

  if (cell_props) {
    intersected_feats <- c("NK","CD8T","CD4T","Bcell","Mono","Neu", intersected_feats)
  }

  if (family) {
    intersected_feats <- c('family_cancer_ind', 'family_cancer_ratio',intersected_feats)
  }

  if (tp53) {
    intersected_feats <- c('dbd','tad','od','nterm','splice','missense','nonsense','frameshift','deletion',intersected_feats)
  }

  return(intersected_feats)

}

##########
# process features
##########
process_features <- function(data,
                             gender,
                             tech,
                             cancer_atdraw,
                             syst_treat,
			     family,
			     tp53) {

  if(gender) {
    data$gender <- ifelse(data$gender == "M", 0, 1)

  }
  if (tech) {
    data$array <- ifelse(data$array == "450", 0, 1)
  }
  if (cancer_atdraw) {
    data$cancer_atdraw <- ifelse(data$cancer_atdraw == 'Yes', 1, 0)
  }

  if (syst_treat) {
    data$systemic_treatment_atdraw <- ifelse(data$systemic_treatment_atdraw == 'Yes', 1, 0)
  }

  if (family) {
    data$family_cancer_ind <- ifelse(data$family_cancer_ind == "cancer", 1, 0)
  }

  if (tp53) {
    lfs_tp53 <- read.csv('Resources/lfs_tp53.csv')
    data$gdna.exon.intron <- lfs_tp53$gdna.exon.intron[match(data$ids,lfs_tp53$ids)]
    data$Variant_Classification <- lfs_tp53$Variant_Classification[match(data$ids,lfs_tp53$ids)]
    data$dbd <- ifelse(data$gdna.exon.intron %in% c("Exon 5","Exon 6","Exon 7","Exon 8","Exon 1-8","Exon 1-11","Exon 7-9"), 1, 0)
    data$tad <- ifelse(data$gdna.exon.intron %in% c("Exon 2","Exon 3","Exon 4", "Exon 1-11", "Exon 1-8"), 1,0)
    data$od <- ifelse(data$gdna.exon.intron %in% c("Exon9","Exon 10","Exon 1-11","Exon 7-8"), 1,0)
    data$nterm <- ifelse(data$gdna.exon.intron %in% c("Exon1","Exon 1-11","Exon 1-8"), 1,0)
    data$missense <- ifelse(data$Variant_Classification == "Missense_Mutation", 1,0)
    data$nonsense <- ifelse(data$Variant_Classification == "Nonsense_Mutation", 1,0)
    data$splice <- ifelse(data$Variant_Classification == "Splice_Site", 1,0)
    data$frameshift <- ifelse(data$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Del_Ins","Frame_Shift_Ins"), 1,0)
    data$deletion <- ifelse(data$Variant_Classification == "Deletion", 1,0)
  }

  return(data)
}

##########
# get outcome
##########
get_outcome <- function(data, age_cutoff) {
  y <- factor(ifelse(data$ageofonset > age_cutoff | is.na(data$ageofonset), "No", "Yes"))
  return(y)
}

##########
# get family variables
##########
# 1) indicator variable if an individaul has at least one family member with cancer
# 2) variable that is the ratio of cancer to healthy individuals
get_family_vars <- function(data) {
  temp <- data %>%      
    group_by(family) %>%
    summarise(num_cancer = sum(cancer_diagnosis !='Unaffected'),
              num_healthy = sum(cancer_diagnosis == 'Unaffected'))
  temp$family_cancer_ind <- ifelse(temp$num_cancer >0, 'cancer', 'no_cancer')
  temp$total_family <- temp$num_cancer + temp$num_healthy
  temp$family_cancer_ratio <- temp$num_cancer/temp$total_family
  temp$num_cancer <- temp$num_healthy <- temp$total_family <- NULL
  data_out <- data.frame(inner_join(temp, data, by='family', multiple = "all"))
  return(data_out)
}


##########b
# predict cancer
##########
# train_dat = train_full
# test_dat = test_full
# age_cutoff = 72
# gender = gender
# tech = tech
# bh_features = remaining_features
pred_cancer_enet <- function(train_dat, 
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			     cancer_atdraw,
			     syst_treat,
			     cell_props,
			     family,
                             tp53,
                             features) {
 
  ind <- grep("array", colnames(test_dat))[1]

  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  train_dat <- process_features(train_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)
  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(train_dat))

  # get y
  train_y <- get_outcome(train_dat, age_cutoff)
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  train_clin <- train_dat[1:ind]
  test_clin <- test_dat[1:ind]
 
  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  N_CV_REPEATS = 2
  nfolds = 5

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
  }
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

pred_cancer_gbm <- function(train_dat,
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
                             cancer_atdraw,
                             syst_treat,
                             cell_props, 
                             family,
                             tp53,
                             features) {

  ind <- grep("array", colnames(test_dat))[1]

  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  train_dat <- process_features(train_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)
  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(train_dat))

  # get y
  train_y <- get_outcome(train_dat, age_cutoff)
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  train_clin <- train_dat[1:ind]
  test_clin <- test_dat[1:ind]

  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]

  NFOLDS <- 5
  folds <- createFolds(train_y, k = NFOLDS, list = TRUE)

  # determines how you train the model.
  NREPEATS <- 2
  fitControl <- trainControl(
    method = "cv",
    number = min(10, NFOLDS),
    classProbs = TRUE,
#    repeats = NREPEATS,
    summaryFunction = twoClassSummary,
#    summaryFunction = auprcSummary,
    allowParallel = TRUE, 
    savePredictions = TRUE,
    index = folds 
  )

  gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 6, 9),
                    n.trees = (0:25)*50, 
                    shrinkage = c(0.001, 0.01, 0.1),
                    n.minobsinnode = 10) 

  gbm.model <- train(x = train_dat
                           , y = train_y
                           , method = "gbm"
#			   , metric = 'AUPRC'
                           , metric = "ROC"
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
# predict cancer
##########
# train_dat = train_full
# test_dat = test_full
# age_cutoff = 72
# gender = gender
# tech = tech
# fam_num = fam_num
# fam_ratio = fam_ratio
# bh_features = remaining_features
pred_cancer_svm <- function(train_dat,
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			     cancer_atdraw,
			     syst_treat,
                             cell_props,
                             family,
                             tp53,
                             features,
			     kernel) {

  ind <- grep("array", colnames(test_dat))[1]

  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  train_dat <- process_features(train_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)
  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(train_dat))

  # get y
  train_y <- get_outcome(train_dat, age_cutoff)
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  train_clin <- train_dat[1:ind]
  test_clin <- test_dat[1:ind]

  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]

  NFOLDS <- 5
  folds <- createFolds(train_y, k = NFOLDS, list = TRUE)

  # determines how you train the model.
  NREPEATS <- 2
  fitControl <- trainControl(
    method = "cv",
    number = min(10, NFOLDS),
    classProbs = TRUE,
#    repeats = NREPEATS,
    summaryFunction = twoClassSummary,
    allowParallel = TRUE,
#    summaryFunction = auprcSummary,
    savePredictions = TRUE,
    index = folds 
  )

  svm.model <- train(x = train_dat
                           , y = train_y
                           , method = kernel
#                           , metric = 'AUPRC'
                           , metric = 'ROC'
                           , trControl = fitControl
                           , verbose = FALSE,
  )
  
  test.predictions <- predict(svm.model 
			      ,newdata = test_dat
			      ,type='prob')

  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))

  return(list(svm.model,temp_dat))

}

##########
# predict cancer
##########
# train_dat = train_full
# test_dat = test_full
# age_cutoff = 72
# gender = gender
# tech = tech
# fam_num = fam_num
# fam_ratio = fam_ratio
# bh_features = remaining_features
pred_cancer_rf <- function(train_dat, 
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			     cancer_atdraw,
			     syst_treat,
                             cell_props,
                             family,
                             tp53,
                             features) {

  ind <- grep("array", colnames(test_dat))[1]

  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  train_dat <- process_features(train_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)
  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(train_dat))

  # get y
  train_y <- get_outcome(train_dat, age_cutoff)
  test_y <- get_outcome(test_dat, age_cutoff)  

  # get clinical data
  train_clin <- train_dat[1:ind]
  test_clin <- test_dat[1:ind]
  
  # get model data
  train_dat <- train_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]

  NFOLDS <- 5
  folds <- createFolds(train_y, k = NFOLDS, list = TRUE)

  # determines how you train the model.
  NREPEATS <- 2
  fitControl <- trainControl(
    method = "cv",
    number = min(10, NFOLDS),
    classProbs = TRUE,
#    repeats = NREPEATS,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary,
#    summaryFunction = auprcSummary,
    savePredictions = TRUE,
    index = folds 
  )

  # mtry: Number of variables randomly sampled as candidates at each split.
  # ntree: Number of trees to grow.
  
  mtry <- sqrt(ncol(train_dat[,colnames(train_dat)]))
  tunegrid <- expand.grid(.mtry=mtry)#, .ntree=c(500,1000,1500,2000),.min.node.size = c(1, 3, 5))
  
  model <- train(x = train_dat
                 , y = train_y
#                 , metric = "AUPRC"
                 , metric = "ROC"
		 , method = "rf"
                 , trControl = fitControl
                 , tuneGrid = tunegrid
                 , importance = TRUE
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

pred_cancer_gbm_test <- function(model,
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			     cancer_atdraw,
			     syst_treat,
                             cell_props,
                             family,
                             tp53,
                             features) {

  ind <- grep("array", colnames(test_dat))[1]
  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(test_dat))

  # get y 
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  test_clin <- test_dat[1:ind]

  # get model data
  test_dat <- test_dat[, intersected_feats]

  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'prob')

  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)

}

pred_cancer_svm_test <- function(model,
                             test_dat,
                             age_cutoff,
                             gender,
                             tech,
			     cancer_atdraw,
			     syst_treat,
                             cell_props,
                             family,
                             tp53,
                             features,
			     kernel) {

  ind <- grep("array", colnames(test_dat))[1]
  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(test_dat))

  # get y
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  test_clin <- test_dat[1:ind]

  # get model data
  test_dat <- test_dat[, intersected_feats]

  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'prob')

  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)

}

pred_cancer_rf_test <- function(model,
			     test_dat,
			     age_cutoff,
                             gender,
                             tech,
			     cancer_atdraw,
			     syst_treat,
                             cell_props,
                             family,
                             tp53,
                             features) {

  ind <- grep("array", colnames(test_dat))[1]
  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(test_dat))

  # get y 
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  test_clin <- test_dat[1:ind]

  # get model data
  test_dat <- test_dat[, intersected_feats]

  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'prob')


  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)

}

pred_cancer_enet_test <- function(model,
			     lambda_index,
                             test_dat,
			     age_cutoff,
                             gender,
                             tech,
			     cancer_atdraw,
			     syst_treat,
                             cell_props,
                             family,
                             tp53,
                             features) {

  ind <- grep("array", colnames(test_dat))[1]
  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(test_dat))

  # get y 
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  test_clin <- test_dat[1:ind]

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

ROCInfo_atopt <- function( data, predict, actual, cost.fp, cost.fn, other_title, model)
{

  data <- data[!duplicated(data$ids),]
  data <- data[!data$tm_donor %in% c("4475"),]

  if (model %in% c("enet","rf","gbm","svmLinear","svmRadial","svmPoly","xgboost","nnet") ) {
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

  if (best_cutoff == "Inf") {
    best_cutoff = 0.5
  }
 
  f1 <- tryCatch(get_f1(data, predict, actual, best_cutoff), error=function(e) return(NA))

  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  auprc <-  PRAUC(y_pred =  data[[predict]], y_true = data[[actual]]) 

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
		auprc = auprc,
                totalcost   = best_cost) )
}

ROCInfo_atcutoff <- function( data, predict, actual, cutoff, other_title, model)
{

  data <- data[!duplicated(data$ids),]
  data <- data[!data$tm_donor %in% c("4475"),]

  if (model %in% c("enet", "rf","gbm","svmLinear","svmRadial","svmPoly","xgboost","nnet") ) {
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
  f1 <- tryCatch(get_f1(data, predict, actual, cutoff), error=function(e) return(NA))

  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  auprc <-  PRAUC(y_pred =  data[[predict]], y_true = data[[actual]])
  
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
                f1 = f1,
		auprc = auprc) )
}

pred_cancer_xgboost <- function(train_dat,
                                test_dat,
                                age_cutoff,
                                gender,
                                tech,
                                cancer_atdraw,
				syst_treat,
                                cell_props,
                                family,
                                tp53,
                                features) {

  ind <- grep("array", colnames(test_dat))[1]
  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  train_dat <- process_features(train_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)
  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(train_dat))

  # get y 
  train_y <- get_outcome(train_dat, age_cutoff)
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  train_clin <- train_dat[1:ind]
  test_clin <- test_dat[1:ind]

  train_dat <- train_dat[intersected_feats]
  test_dat <- test_dat[intersected_feats]

  NFOLDS <- 5
  folds <- createFolds(train_y, k = NFOLDS, list = TRUE)

  # determines how you train the model.
  NREPEATS <- 2
  fitControl <- trainControl(
    method = "cv",
    number = min(10, NFOLDS),
    classProbs = TRUE,
#    repeats = NREPEATS,
    allowParallel = TRUE,
#    summaryFunction = auprcSummary,
    summaryFunction = twoClassSummary,
    savePredictions = TRUE,
    index = folds 
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
#                     , metric = "AUPRC"
                     , metric = "ROC"
                     , tuneGrid = gbmGrid
                     , trControl = fitControl
                     , verbose = FALSE
  )

  test.predictions <- predict(gbm.model
                              ,newdata = test_dat
                              ,type='prob')

  print(test.predictions)

  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))

  return(list(gbm.model,temp_dat))

}

pred_cancer_nnet <- function(train_dat,
                                test_dat,
                                age_cutoff,
                                gender,
                                tech,
                                cancer_atdraw,
                                syst_treat,
                                cell_props,
				family,
				tp53,
                                features) {

  ind <- grep("array", colnames(test_dat))[1]
  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  train_dat <- process_features(train_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)
  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(train_dat))
  # get y 
  train_y <- get_outcome(train_dat, age_cutoff)
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  train_clin <- train_dat[1:ind]
  test_clin <- test_dat[1:ind]

  train_dat <- train_dat[intersected_feats]
  test_dat <- test_dat[intersected_feats]

  NFOLDS <- 5
  folds <- createFolds(train_y, k = NFOLDS, list = TRUE)

  # determines how you train the model.
  NREPEATS <- 2
  fitControl <- trainControl(
    method = "cv",
    number = min(10, NFOLDS),
    classProbs = TRUE,
#    repeats = NREPEATS,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary,
#    summaryFunction = auprcSummary,
    savePredictions = TRUE,
    index = folds 
  )

  nnetGrid <- expand.grid(
	.decay = c(0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5), 
	.size = c(3, 5, 10, 20)
  )

  nnet.model <- train(x = train_dat
                     , y = train_y
                     , method = "nnet"
                     , tuneGrid = nnetGrid
                     , trControl = fitControl
                     , verbose = TRUE
#		      , metric = "AUPRC"
                     , metric = "ROC"
		     , maxit = 100
		     )

  test.predictions <- predict(nnet.model
                              ,newdata = test_dat
                              ,type='prob')

  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))

  return(list(nnet.model,temp_dat))
   
}

pred_cancer_nnet_test <- function(model,
                                  test_dat,
                                  age_cutoff,
                                  gender,
                                  tech,
                                  cancer_atdraw,
                                  syst_treat,
                                  cell_props,
				  family,
                                  tp53,
                                  features) {

  ind <- grep("array", colnames(test_dat))[1]
  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(test_dat))

  # get y 
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  test_clin <- test_dat[1:ind]

  # get model data
  test_dat <- test_dat[, intersected_feats]

  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'prob')

  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)

}

pred_cancer_xgboost_test <- function(model,
                                         test_dat,
                                         age_cutoff,
                                         gender,
                                         tech,
                                         cancer_atdraw,
					 syst_treat,
	                                 cell_props,
					 family,
                                         tp53,
	                                 features) {

  ind <- grep("array", colnames(test_dat))[1]
  intersected_feats <- get_features(gender, tech, cancer_atdraw, syst_treat,cell_props, family, tp53, features)

  test_dat <- process_features(test_dat, gender, tech, cancer_atdraw, syst_treat, family, tp53)

  intersected_feats <- intersect(intersected_feats,colnames(test_dat))

  # get y 
  test_y <- get_outcome(test_dat, age_cutoff)

  # get clinical data
  test_clin <- test_dat[1:ind]

  # get model data
  test_dat <- test_dat[, intersected_feats]

  test.predictions <- predict(model,
                              data.matrix(test_dat),
                              type = 'prob')

  # combine predictions and real labels
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions, test_label = test_y, test_clin))
  return(temp_dat)

}

platt_scaling <- function(train_results, test_results,model) {
  actual <- "test_label"
  if (model %in% c("svmPoly", "rf","svmLinear","xgboost","svmRadial","gbm","nnet") ) {
    predict <- "test_pred.Yes"
  } else {
    predict <- "test_pred"
  }
  #calculating Log Loss without Platt Scaling
  ll <- LogLoss(train_results[[predict]],train_results[[actual]])
  ll_df <- data.frame(x=train_results[[predict]],y=as.factor(train_results[[actual]]))
  model_log <-glm(y~x,data = ll_df,family = binomial)
  saveRDS(model_log,'rds/recalibration_model.rds')
  #predicting on the cross validation after platt scaling
  result_platt<-predict(model_log,ll_df["x"],type = "response")
  train_results$test_pred_calibrated.Yes <- result_platt
  
  # Predicting on the test dataset using Platt Scaling
  ll_df_test<-data.frame(x=test_results[[predict]])
  result_test_platt<-predict(model_log,ll_df_test,type="response")
  test_results$test_pred_calibrated.Yes <- result_test_platt
  
  return(list(train_results,test_results))
}

get_all_auc <- function(dir,ext) {
  auc <- list() ; sens <- list() ; spec <- list() ; f1 <- list() ; auprc <- list()
  all_models <- list.files(dir,ext)
  for (model in all_models) {
    cat(paste0(model,'\n'))
    r_all <- readRDS(paste0(dir,model))
    id <- gsub("_ageofonset_ROCInfoTest.rds|_ageofonset_ROCInfoVal.rds","",model)
    auc[id] <- r_all[[6]]
    sens[id] <- r_all[[7]]
    spec[id] <- r_all[[8]]
    f1[id] <- r_all[[9]]
    auprc[id] <- r_all[[10]]

  }
  auc_all <- data.frame(auc = do.call(rbind,auc),
                         sens = do.call(rbind,sens),
                         spec = do.call(rbind,spec),
                         f1 = do.call(rbind,f1),
			 auprc = do.call(rbind,auprc))
  auc_all$preprocessing <- rownames(auc_all)
  return(auc_all)
}


get_avg_results <- function(dir) {
  files <- list.files(dir, pattern = "_auc_all.csv", recursive = TRUE, full.names = TRUE)
  all <- do.call(rbind, lapply(files, read.csv, row.names = 'X'))
  ag <- aggregate(. ~ dataset + preprocessing, all, function(x) c(mean = mean(x,na.rm=TRUE), sd = sd(x,na.rm=TRUE)))
  write.csv(ag,file.path(dir,'all_avg_model_results.csv'),row.names=F, quote=F)
}


