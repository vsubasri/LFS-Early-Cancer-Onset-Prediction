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
parser$add_argument("-v","--value", action="store")
parser$add_argument("-s", "--sex", action="store_true", default=FALSE, help="Remove probes in sex chromosomes")
value <- args$value

#####################################################################
############################ Main Script ############################
#####################################################################

cat("[ Reading in input sample sheets ]","\n")

#cases batch1
id_map_tor <- read.csv('Data/LFS_Mut/450k/cases_toronto/SampleSheet.csv', stringsAsFactors = F, fileEncoding="UTF-8-BOM")

#cases batch2
id_map_mon <- read.csv('Data/LFS_Mut/450k/cases_montreal/SampleSheet.csv', stringsAsFactors = F, fileEncoding="UTF-8-BOM")

# combine id_map and id_map_other
id_tor_mon <- rbind(id_map_tor, id_map_mon) ; rm(id_map_tor, id_map_mon)

id_cases <- data.frame(Project=id_tor_mon$Project,
		       SentrixID=paste0(id_tor_mon$Sentrix_ID,"_",id_tor_mon$Sentrix_Position), 
                       ids=cbind(lapply(strsplit(id_tor_mon$Sample_Name,'#'), function(x) x[length(x)])))

id_map_utah <- read.csv('Data/LFS_Mut/850k/schiffman/SUB12648-71.csv', stringsAsFactors = F)
id_utah <- data.frame(Project=id_map_utah$Project,
		      SentrixID=paste0(id_map_utah$Sentrix_ID,"_",id_map_utah$Sentrix_Position),
                      ids=id_map_utah$Sample_Name)

id_map_batch5 <- read.csv('Data/LFS_Mut/850k/batch5/SUB14749-51.csv', stringsAsFactors = F)
id_batch5 <- data.frame(Project=id_map_batch5$Project,
                      SentrixID=paste0(id_map_batch5$Sentrix_ID,"_",id_map_batch5$Sentrix_Position),
                      ids=id_map_batch5$Sample_Name)


id_map <- rbind(id_cases, id_utah, id_batch5) %>% getIdName() 

clin <- read.csv('Data/clin_data/lfs_mut_clinical_comprehensive.csv', stringsAsFactors = F, fileEncoding="UTF-8-BOM")
id_map_clin <- merge(id_map,clin, by.x="ids", by.y="sample",all.x=T) 

###################################################################################################################
# Preprocess samples
###################################################################################################################

cat("[ Process 450k ]","\n")
data_450k_raw <- processData('Data/LFS_Mut/450k',id_map_clin,"raw",value,"450")
data_450k_noob <- processData('Data/LFS_Mut/450k',id_map_clin,"noob",value,"450")
saveRDS(data_450k_noob,paste0('rds/450k_Noob_',value,'.rds'))

cat("[ Process 850k Schiffman ]","\n")
data_schiffman_raw <- processData('Data/LFS_Mut/850k/schiffman',id_map_clin,"raw",value,"850")
data_schiffman_noob <- processData('Data/LFS_Mut/850k/schiffman',id_map_clin,"noob",value,"850")
saveRDS(data_schiffman_noob,paste0('rds/Schiffman_Noob_',value,'.rds'))

cat("[ Process 850k batch5]","\n")
data_batch5_raw <- processData('Data/LFS_Mut/850k/batch5',id_map_clin,"raw",value,"850")
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

