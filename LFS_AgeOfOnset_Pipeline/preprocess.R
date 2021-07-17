#####################################################################
############################# Libraries #############################
#####################################################################

require(bumphunter)
require(minfi)
require(dplyr)
require(IlluminaHumanMethylation450kmanifest)
require(IlluminaHumanMethylationEPICmanifest)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
require(ggplot2)
require(reshape2)
require(argparse)
require(ggplot2)
require(reshape2)

setwd('/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/')
source('Scripts/util_functions.R')

parser <- ArgumentParser()
parser$add_argument("-s","--samplesheet", action="store")
parser$add_argument("-d","--dir", action="store")
samplesheet <- args$samplesheet
dir <- args$dir

#####################################################################
############################ Main Script ############################
#####################################################################

cat("[ Reading in input sample sheets ]","\n")

id_map_clin <- read.csv(samplesheet, stringsAsFactors = F, fileEncoding="UTF-8-BOM")
id_map_450k <- 

###################################################################################################################
# Preprocess samples
###################################################################################################################

cat("[ Process 450k ]","\n")
data_450k_noob <- processData(dir,id_map_clin,"noob","beta","450")
saveRDS(data_450k_noob,paste0('rds/Noob_',value,'.rds'))

cat("[ Process 850k Schiffman ]","\n")
data_batch5_noob <- processData('Data/LFS_Mut/850k/batch5',id_map_clin,"noob",value,"850")
saveRDS(data_batch5_noob,paste0('rds/Batch5k_Noob_',value,'.rds'))

cols <- Reduce(intersect, list(names(data_450k_noob),names(data_schiffman_noob),names(data_batch5_noob)))
data_450k_raw <- data_450k_raw[cols] ; data_schiffman_raw <- data_schiffman_raw[cols] ; data_batch5_raw <- data_batch5_raw[cols]
data_450k_noob <- data_450k_noob[cols] ; data_schiffman_noob <- data_schiffman_noob[cols] ; data_batch5_noob <- data_batch5_noob[cols]

cat("[ Combine data ]","\n")
dataRaw <- rbind(data_450k_raw,data_schiffman_raw,data_batch5_raw) ; saveRDS(dataRaw,paste0('rds/Raw_',value,'.rds'))
dataNoob <- rbind(data_450k_noob,data_schiffman_noob,data_batch5_noob) ; saveRDS(dataNoob,paste0('rds/Noob_',value,'.rds'))

Raw_cleaned <- dataRaw[45:length(dataRaw)]
Raw_cleaned <- Raw_cleaned[, colSums(is.na(Raw_cleaned)) == 0]
Noob_cleaned <- dataNoob[45:length(dataNoob)]
Noob_cleaned <- Noob_cleaned[, colSums(is.na(Noob_cleaned)) == 0]
pheno <- dataNoob[1:44]

########################## PLOT RAW DATA ###########################

indices <- which(is.na(pheno$cancer_diagnosis))

pc <- prcomp(as.matrix(Raw_cleaned[-indices,]), scale = TRUE)
pc_clin <- cbind(pheno[-indices,],pc$x)
write.csv(pc_clin,paste0('Output/Raw_',value,'_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0('Raw_',value,'_PCA_summary.csv'))
generate_pcplots(pc_clin,paste0('Raw_',value))

########################## PLOT NOOB DATA ##########################

pc <- prcomp(as.matrix(Noob_cleaned[-indices,]), scale = TRUE)
pc_clin <- cbind(pheno[-indices,],pc$x)
write.csv(pc_clin,paste0('Output/Noob_',value,'_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0('Noob_',value,'_PCA_summary.csv')
generate_pcplots(pc_clin,paste0('Noob_',value))

cat("[ Plot technical replicates ]","\n")

duplicated_data <- get_technicalreplicates(dataNoob)
plot_concordance(duplicated_data,paste0('Noob_',value))

pc <- data.frame(duplicated_data[45:length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:44],pc$x)
write.csv(pc_clin,paste0('Output/Noob_',value,'_TechnicalReplicates_PCA.csv',quote=F,row.names=F))

generate_pcsummary(pc_clin,paste0('Noob_',value,'_TechnicalReplicates_PCA_summary.csv'))
generate_pcplots(pc_clin,paste0('Noob_',value,'_TechnicalReplicates'))

########################## REMOVE SEX CHR PROBES ##########################

if (args$sex) {

    Noob_nosex <- remove_sex(Noob_cleaned, pheno, paste0('Noob_',value)
 
    pc <- prcomp(as.matrix(Noob_nosex[-indices,]), scale = TRUE)
    pc_clin <- cbind(pheno[-indices,],pc$x)
    write.csv(pc_clin,paste0('Output/Noob_',value,'_NoSex_PCA.csv'),quote=F,row.names=F)

    generate_pcsummary(pc_clin,paste0('Noob_',value,'_NoSex_PCA_summary.csv')
    generate_pcplots(pc_clin,paste0('Noob_',value,'_NoSex'))

}

