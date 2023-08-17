## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Valli Subasri
##
## Date Created: 2019-04-11
##
## Email: vallisubasri@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## load up the packages ## 
suppressMessages(require(ggplot2))
suppressMessages(require(argparse))
suppressMessages(require(dplyr))
suppressMessages(require(limma))
suppressMessages(require(data.table))
suppressMessages(require(umap))
suppressMessages(require(tidyr))

## initiate parser ##
parser <- ArgumentParser()
parser$add_argument("--datafile", action="store")
parser$add_argument("--value", action="store")
parser$add_argument("--covar", action="store")
parser$add_argument("--seed", action="store")
parser$add_argument("--nsplit", action="store",type="integer")
parser$add_argument("--outdir", action="store")

## set variables from parser
args <- parser$parse_args()
datafile <- args$datafile
value <- args$value
covar <- args$covar
seed <- args$seed
nsplit <- args$nsplit
outdir <- args$outdir

## set seed
set.seed(seed)
## set working directory and utility functions ##
setwd('/hpf/largeprojects/davidm/vsubasri/methyl_data')
source('Scripts/LFS_ageofonset/util_functions.R')

#####################################################################
############################# Functions #############################
#####################################################################

## function to identify 850k replicate with the closest age of sample collection ##
get_850replicate_indices <- function(data_450k,data_850k) {
  age_450 <- as.numeric(data_450k["agesamplecollection"])
  replicate_850 <- data_850k[data_850k$tm_donor == data_450k["tm_donor"],]
  replicate_850$agediff <- replicate_850$agesamplecollection - age_450
  min <- replicate_850$ids[which.min(replicate_850$agediff)]
  return(min)
}

## identify match technical replicates between 450k and 850k ##
get_matching_replicates <- function(data) {
  cat("[ Getting technical replicates ]","\n")
  data <- data[!is.na(data$tm_donor),]
  data_450k <- data[data$array == "450",]
  data_850k <- data[data$array == "850",]
  tech_450k <- data_450k[data_450k$tm_donor %in% data_850k$tm_donor,]
  tech_450k <- tech_450k[!duplicated(tech_450k$tm_donor),]
  tech_ids <- apply(tech_450k,1, get_850replicate_indices,data_850k)
  tech_850k <- data_850k[data_850k$ids %in% tech_ids,]
  tech_450k <- tech_450k[tech_450k$tm_donor %in% tech_850k$tm_donor,]
  tech_450k <- tech_450k[match(tech_850k$tm_donor,tech_450k$tm_donor),]
  tech_all <- list(tech_450k,tech_850k)
  return(tech_all)
}


## project 450k data onto the 850k space ##
project_450k <- function(data,tech_all,covar,outdir,value,type=NULL) {
  cat("[ Project 450k data onto 850k space ]","\n")
  data_450k <- data[!is.na(data$cancer_diagnosis),]
  tech_450k <- tech_all[[1]]
  tech_850k <- tech_all[[2]]

  ind <- grep("cg", colnames(data))[1]-1

  models <- list()
  if (!is.null(type)) {
    data_450k <- data_450k[data_450k$cancerstatus == type,]
    tech_450k <- tech_450k[tech_450k$cancerstatus == type,]
    tech_850k <- tech_850k[tech_850k$cancerstatus == type,]
  } 

  cols <- colnames(tech_450k)[(ind+1):length(tech_450k)]; lm_all <- list()
  for (i in 1:length(cols)) {
    probe <- cols[i]
    probe_450=tech_450k[,probe]
    probe_850=tech_850k[,probe]
    age_sample_collection=tech_450k$agesamplecollection
    gender=tech_450k$gender
    new_data <- data.frame(probe_450 = data_450k[,probe],
			   age_sample_collection = data_450k$agesamplecollection,
			   gender = data_450k$gender)
    if (covar) {
	models[[probe]] <- lm(probe_850 ~ probe_450 + age_sample_collection + gender)
	lm_probe <- predict(lm(probe_850 ~ probe_450 + age_sample_collection + gender),new_data)
    } else {
	lm_probe <- predict(lm(probe_850 ~ probe_450),new_data)
    } 
    lm_all[[probe]] <-lm_probe
    if (i %% 50000 == 0) {
      print(i)
    }
  }
  saveRDS(models,paste0(outdir,'rds/arrayEffectRemoval_',value,'.rds'))
  corrected_meth <- do.call(cbind,lm_all)
  corrected_data <- cbind(data_450k[1:ind],corrected_meth)
  return(corrected_data)
}


## run step2 of correction (legacy function) ##
run_step2 <- function(data_450k,data_850k,corrected_450k) {
  ## Step 2a) predict cancer signal in corrected 450k
  step2_450k <- model_cancer_signal2(data_450k,corrected_450k)
  saveRDS(step2_450k,'rds/NoobFinal450k2_beta_NoSex.rds')
  pc <- prcomp(~ ., data=as.matrix(step2_450k[(ind+1):length(step2_450k)]), na.action=na.omit, scale = TRUE)
  pc_clin <- cbind(step2_450k[1:ind],pc$x)
  write.csv(pc_clin,paste0('Output/NoobFinal450k2_beta_NoSex_PCA.csv'),quote=F,row.names=F)
  generate_pcsummary(pc_clin,"NoobFinal450k2_beta_NoSex_PCA_summary.csv")
  generate_pcplots(pc_clin,"NoobFinal450k2_beta_NoSex")

  ## Step 2b) predict cancer signal in 850k 
  step2_850k <- model_cancer_signal(data_450k,data_850k)
  saveRDS(step2_850k,'rds/NoobFinal850k_beta_noSex.rds')
  pc <- prcomp(as.matrix(step2_850k[(ind+1):length(step2_850k)]), scale = TRUE)
  pc_clin <- cbind(step2_850k[1:ind],pc$x)
  write.csv(pc_clin,paste0('Output/NoobFinal850k_beta_NoSex_PCA.csv'),quote=F,row.names=F)
  generate_pcsummary(pc_clin,"NoobFinal850k_beta_NoSex_PCA_summary.csv")
  generate_pcplots(pc_clin,"NoobFinal850k_beta_NoSex")

  step2_data <- rbind(step2_450k,step2_850k)
  saveRDS(step2_data,'rds/NoobFinal_beta_NoSex.rds')
  return(step2_data)
}

## run step 1 of correction with cancer/unaffected separately ##
run_step1 <- function(data,value,covar,outdir) {
  tech_all <- get_matching_replicates(data)
  data_450k <- data[data$array == "450",]
  u <- project_450k(data_450k,tech_all,covar,outdir,value,"Unaffected")
  saveRDS(u,paste0(outdir,"rds/NoobCorrected_Null_",value,"_NoSex.rds"))
  c <- project_450k(data_450k,tech_all,outdir,value,"Cancer")
  saveRDS(c,paste0(outdir,"rds/NoobCorrected_Cancer_",value,"_NoSex.rds"))
  u <- u[colnames(u)[colnames(u) %in% colnames(c)]]
  c <- c[colnames(c)[colnames(c) %in% colnames(u)]]
  corrected_data <- rbind(u,c)
  saveRDS(corrected_data,paste0(outdir,"rds/NoobCorrected450k_",value,".rds"))
  return(corrected_data)
}

## run step 1 of correction together ##
run_step1together <- function(data,value,covar,outdir) {
  tech_all <- get_matching_replicates(data)
  data_450k <- data[data$array == "450",]
  c <- project_450k(data_450k,tech_all,covar,outdir,value)
  saveRDS(c,paste0(outdir,"rds/NoobCorrected450k_",value,".rds"))
  return(c)
}

#####################################################################
############################ Main Script ############################
#####################################################################

## read in input data ##
cat("[ Reading in input data ]","\n")
data <- readRDS(datafile)
data <- data[!is.na(data$cancer_diagnosis),]
data$array <- factor(data$array)

ind <- grep("cg", colnames(data))[1]-1

## remove outliers ## 
cat("[ PCA of all data before correction ]","\n")
pc <- prcomp(as.matrix(data[(ind+1):length(data)]), scale = TRUE)
pc_clin <- cbind(data[1:ind],pc$x)
keep <- remove_outliers(pc_clin,3)
data <- data[data$SentrixID %in% keep,]
data_450k <- data[data$array == "450",]

## correct 450k data ##
corrected_450k <- run_step1together(data,value,covar,outdir)
data_850k <- data[data$array == "850",]
all_corrected <- rbind(corrected_450k,data_850k)
all_corrected <- all_corrected[!is.na(all_corrected[length(all_corrected)]),]
saveRDS(all_corrected,paste0(outdir,"rds/NoobCorrected_",value,".rds"))

## plot correction ##
cat("[ PCA of all data after correction ]","\n")
pc <- prcomp(as.matrix(all_corrected[(ind+1):length(all_corrected)]), scale = TRUE)
pc_clin <- cbind(all_corrected[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/NoobCorrected_',value,'_PCA.csv'),quote=F,row.names=F)
generate_pcsummary(pc_clin,paste0("NoobCorrected_",value,"_PCA_summary.csv"),outdir)
generate_pcplots(pc_clin,paste0("NoobCorrected_",value),outdir)

u <- umap(all_corrected[(ind+1):length(all_corrected)]) ; ud <- cbind(all_corrected[1:ind],u$layout)
write.csv(ud,paste0(outdir,'Output/Umap_NoobCorrected_',value,'.csv'))

## plot correction for technical replicates only ##
cat("[ PCA of all technical replicates after correction ]","\n")
duplicated_data <- get_technicalreplicates(all_corrected)
plot_concordance(duplicated_data,paste0("NoobCorrected_",value),outdir)

pc_beta <- data.frame(duplicated_data[(ind+1):length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:ind],pc$x)
write.csv(pc_clin,paste0(outdir,'Output/NoobCorrected_',value,'_TechnicalReplicates_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0('NoobCorrected_',value,'_TechnicalReplicates_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0('NoobCorrected_',value,'_TechnicalReplicates'),outdir)

## generate train/test/val splits
datasplits <- getTrainTestVal(all_corrected,seed,nsplit)

data_train <- datasplits[[1]]
saveRDS(data_train,paste0(outdir,"rds/NoobCorrected_",value,"_TrainingSet.rds"))
data_test <- datasplits[[2]]
saveRDS(data_test,paste0(outdir,"rds/NoobCorrected_",value,"_TestSet.rds"))
data_val <- datasplits[[3]]
saveRDS(data_val,paste0(outdir,"rds/NoobCorrected_",value,"_ValidationSet.rds"))


