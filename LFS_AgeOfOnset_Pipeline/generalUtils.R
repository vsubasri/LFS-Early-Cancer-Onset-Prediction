#####################################################################
############################# Libraries #############################
#####################################################################

require(dplyr)
require(reshape2)
require(minfi)

#####################################################################
############################# Functions #############################
#####################################################################

##########
# Function to preprocess raw methylation data
##########

preprocessMethod <- function(data, preprocess) {
  if (preprocess == 'raw') {
    Mset <- preprocessRaw(data)
  }
  if (preprocess == 'quan') {
    Mset <- preprocessQuantile(data, fixOutliers = TRUE,
                               removeBadSamples = TRUE, badSampleCutoff = 10.5,
                               quantileNormalize = TRUE, stratified = TRUE,
                               mergeManifest = FALSE, sex = NULL)
  }
  if (preprocess == 'illumina') {
    Mset  <- preprocessIllumina(data)
  }
  if (preprocess == 'swan') {
    Mset  <-preprocessSWAN(data)
  }
  if (preprocess == 'funnorm') {
    Mset <-preprocessFunnorm(data)
  }
  if (preprocess == 'noob') {
    Mset <-preprocessNoob(data, dyeMethod="single")
  }
  return(Mset)
}

##########
# Function to remove batch effect
##########

remove_batch_effect <- function(data, outdir,id) {
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
# Function that preprocesses data by first converts 850 arrays to to 450 arrays, preprocesssimputes sex 
##########

processData <- function(directory,id_map_clin,method,valtype,arraytype) {
	rgCases <- read.metharray.exp(directory, recursive = T, force=TRUE)
	rgCases  <- convertArray(rgCases,outType = c("IlluminaHumanMethylation450k"))
	preprocessedCases <- preprocessMethod(rgCases,method)
	pheno <- id_map_clin[match(rownames(preprocessedCases@colData),id_map_clin$SentrixID),]
	pheno$cancerstatus <- ifelse(pheno$cancer_diagnosis != "Unaffected","Cancer","Unaffected")
	pheno$array <- arraytype
	GRset <- mapToGenome(preprocessedCases)
	estSex <- getSex(GRset)
	if (valtype == "beta") {
		meth <- data.frame(getBeta(GRset))
	} else if (valtype == "M") {
		meth <- data.frame(getM(GRset))
	} else {
		cat("Invalid methylation type value")
	}
	pheno$gender <- estSex[, "predictedSex"]
	data <- cbind(pheno,t(meth))
	return(data)
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
# Function that gets malkin id from sentrix id
##########

getIdName <- function(data) {
  column_split <- strsplit(as.character(data$ids), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  sub_ids <- gsub('-', '', sub_ids)
  data$ids <- sub_ids
  data$identifier <- NULL
  return(data)
}


##########
# Function that splits data into training, test and validation set
##########

getTrainTestVal <- function(data,nseed,nsamps) {
	set.seed(nseed)
	ind <- 44
  excl_ids <- c(2765,2957,3365,3367,3503,3885,4033,3301,3692,3733,3884,3886,3887,3888,3889,3892,4036,4037,4164,4165,4169,3304)
	data <- data[!data$ids %in% excl_ids,]
	## create test set
  schiff_ids <- seq(6056,6147,1)
	data_val <- data[data$ids %in% schiff_ids,]
	## create validation set
	data <- data[!data$ids %in% schiff_ids,]
	test_ids <- data[1:ind] %>%
        	filter(!ids %in% schiff_ids) %>%
        	group_by(ageoutcome = ifelse(data$ageofonset > 12*6 | is.na(data$ageofonset),"GT6","LT6")) %>%
        	sample_n(nsamps) %>%
        	select(tm_donor)
	data_test <- data[data$tm_donor %in% test_ids$tm_donor,]
	## create training set
	data_train <- data[!data$tm_donor %in% c(data_test$tm_donor,data_val$tm_donor),]

	## remove outliers in training data
	p_train <- prcomp(data_train[(ind+1):length(data_train)],scale=T)
	p_train <- cbind(data_train[1:ind],p_train$x)
	keep_train <- remove_outliers(p_train,3)
	data_train <- data_train[data_train$SentrixID %in% keep_train,]
	return(list(data_train,data_test,data_val))
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
# Function to remove loci overlapping SNPs with MAF of choice and cross reactive SNPs
##########

remove_SNPs_xRP <- function(GRset,pheno,output) {
  GRset_noSNPs <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)
  B_noSNPs <- data.frame(getBeta(GRset_noSNPs))
  B_noSNPs <- B_noSNPs[rowSums(is.na(B_noSNPs)) == 0,]
  data_B_noSNPs <- cbind(pheno,t(B_noSNPs))
  saveRDS(data_B_noSNPs,paste0("rds/",output,"_noSNPs.rds"))

  xReactiveProbes <- read.csv(file="Resources/cross_reactive_probes.csv", stringsAsFactors=FALSE)
  keep <- !(featureNames(GRset_noSNPs) %in% xReactiveProbes$TargetID)
  GRset_noSNPs_noxRP <- GRset_noSNPs[keep,]
  B_noSNPs_noxRP <- data.frame(getBeta(GRset_noSNPs_noxRP))
  B_noSNPs_noxRP <- B_noSNPs_noxRP[rowSums(is.na(B_noSNPs_noxRP)) == 0,]
  data_B_noSNPs_noxRP <- cbind(pheno,t(B_noSNPs_noxRP))
  saveRDS(data_B_noSNPs_noxRP,paste0("rds/",output,"_noSNPs_noxRP.rds"))
}


##########
# Function to remove probes located in sex chromosomes
##########

remove_sex <- function(cleaned,pheno,output) {
  sexprobes <- scan('Scripts/sexprobes.txt',what="character")
  nosex <- cleaned[, !colnames(cleaned) %in% sexprobes]
  data_nosex <- cbind(pheno,nosex)
  saveRDS(data_nosex,paste0("rds/",output,"_NoSex.rds"))
  return(nosex)
}

##########
# Function that gets all samples profiled on the 450k and 850k array
##########

get_technicalreplicates <- function(data) {
  duplicated_ids_diffarray = data[c("tm_donor","array")] %>%
    filter(!is.na(tm_donor)) %>%
    group_by(tm_donor) %>%
    mutate(types = n_distinct(array)) %>%
    filter(types > 1) %>%
    distinct(tm_donor, .keep_all = TRUE) %>%
    pull(tm_donor)
  duplicated_data <- data[data$tm_donor %in% duplicated_ids_diffarray,] %>%
    distinct(tm_donor, array, .keep_all = TRUE)
  return(duplicated_data)
}

##########
# Function that samples 100000 probes and plot the difference in methylation alue across technical replicates
##########

plot_concordance <- function(duplicated_data,value,outdir) {
  allprobes = colnames(duplicated_data)[grepl("cg",colnames(duplicated_data))]
  sample_probes = sample(allprobes,100000,replace=FALSE)
  sample_beta = duplicated_data[c("tm_donor","array",sample_probes)]
  plot_dups <- melt(sample_beta,measure.vars=sample_probes, id.vars= c("tm_donor","array")) %>%
  dcast(tm_donor + variable ~ array)
  tech_corr <- round(cor(plot_dups$`450`,plot_dups$`850`,method="spearman"),digits=3)
  plot_dups$diff <- abs(plot_dups$`450`-plot_dups$`850`)
  distsum <- sum(plot_dups[, 'diff'])

  pdf(paste0(outdir,'Plots/',value,'_450_850_concordance.pdf'),width=9,height=7)
  cplot <- ggplot(plot_dups,aes(x=`450`,y=`850`)) +
    geom_point(aes(colour = diff)) +
    theme_bw() +
    ggtitle(paste0("Correlation = ",tech_corr,", Distance = ",distsum)) +
    scale_colour_gradient2(low="yellow",mid="black",high="blue") +
    labs(x="450k Array",y="850k Array",color=value)
  print(cplot)
  suppressMessages(dev.off())
}

##########
# Function to identify 850k replicate with the closest age of sample collection 
##########

get_850replicate_indices <- function(data_450k,data_850k) {
  age_450 <- as.numeric(data_450k["agesamplecollection"])
  replicate_850 <- data_850k[data_850k$tm_donor == data_450k["tm_donor"],]
  replicate_850$agediff <- replicate_850$agesamplecollection - age_450
  min <- replicate_850$ids[which.min(replicate_850$agediff)]
  return(min)
}

##########
# Function to identify match technical replicates between 450k and 850k
##########

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

##########
# Function to project 450k data onto the 850k space
##########

project_450k <- function(data,tech_all,covar,type=NULL) {
  cat("[ Project 450k data onto 850k space ]","\n")
  data_450k <- data[!is.na(data$cancer_diagnosis),]
  tech_450k <- tech_all[[1]]
  tech_850k <- tech_all[[2]]

  if (!is.null(type)) {
    data_450k <- data_450k[data_450k$cancerstatus == type,]
    tech_450k <- tech_450k[tech_450k$cancerstatus == type,]
    tech_850k <- tech_850k[tech_850k$cancerstatus == type,]
  } 

  cols <- colnames(tech_450k)[45:length(tech_450k)]; lm_all <- list()
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
  lm_probe <- predict(lm(probe_850 ~ probe_450 + age_sample_collection + gender),new_data)
    } else {
  lm_probe <- predict(lm(probe_850 ~ probe_450),new_data)
    } 
    lm_all[[probe]] <-lm_probe
    if (i %% 50000 == 0) {
      print(i)
    }
  }
  corrected_meth <- do.call(cbind,lm_all)
  corrected_data <- cbind(data_450k[1:44],corrected_meth)
  return(corrected_data)
}


##########
# Function to run correction
##########

run_correction <- function(data,value,covar,outdir) {
  tech_all <- get_matching_replicates(data)
  data_450k <- data[data$array == "450",]
  c <- project_450k(data_450k,tech_all,covar)
  saveRDS(c,paste0(outdir,"rds/NoobCorrected450k_",value,".rds"))
  return(c)
}

