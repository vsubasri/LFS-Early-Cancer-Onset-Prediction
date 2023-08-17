library(bumphunter)
library(dplyr)
library(minfi)
library(doParallel)

registerDoParallel(cores = 3)

DELTA_BETA_THRESH = .10 # DNAm difference threshold
NUM_BOOTSTRAPS = 10   # number of randomizations

ind <- 67

#cancertypes: string with cancer types of interest each separated by | with no spaces before or after
filter_data <- function(data, age_upper, age_lower, cancertypes) {
  types <- unlist(strsplit(cancertypes, "|", fixed=TRUE))
  for (i in 1:length(types)) {
    data <- data[data$cancer_diagnosis_diagnoses %in% types,]
  }
  if (!missing(age_upper) & !missing(age_lower)) {
    data <- subset(data, data$age_sample_collection < age_upper & 
                   data$age_sample_collection > age_lower)
  }
  return(data)
}

#preprocess data to include samples of interest prior to calling this function
process_params <- function(data, ratio_set) {
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

bump_hunter <- function(params, designMatrix, name) {
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
  
  tab <- bumphunter(params[[2]], 
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
  
  save(tab, file = paste0('rds/',name,".rds"), compress = "xz")
  
}
