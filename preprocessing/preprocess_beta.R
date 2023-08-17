suppressMessages(require(bumphunter))
suppressMessages(require(minfi))
suppressMessages(require(dplyr))
suppressMessages(require(IlluminaHumanMethylation450kmanifest))
suppressMessages(require(IlluminaHumanMethylationEPICmanifest))
suppressMessages(require(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
suppressMessages(require(ggplot2))
suppressMessages(require(reshape2))

setwd('/hpf/largeprojects/davidm/vsubasri/methyl_data')
source('Scripts/LFS_ageofonset/util_functions.R')

################################ READ IN INPUT ID FILES #################################

outdir <- '/hpf/largeprojects/davidm/vsubasri/methyl_data/'

cat("[ Reading in input sample sheets ]","\n")

#cases batch1
id_map_tor <- read.csv('Data/LFS_Mut/450k/cases_toronto/SampleSheet.csv', stringsAsFactors = F, fileEncoding="UTF-8-BOM")

#cases batch2
id_map_mon <- read.csv('Data/LFS_Mut/450k/cases_montreal/SampleSheet.csv', stringsAsFactors = F, fileEncoding="UTF-8-BOM")

#combine id_map and id_map_other
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
data_450k_raw <- processData('Data/LFS_Mut/450k',id_map_clin,"raw","beta","450")
data_450k_noob <- processData('Data/LFS_Mut/450k',id_map_clin,"noob","beta","450")
saveRDS(data_450k_noob,'rds/450k_Noob_beta.rds')

cat("[ Process 850k Schiffman ]","\n")
data_schiffman_raw <- processData('Data/LFS_Mut/850k/schiffman',id_map_clin,"raw","beta","850")
data_schiffman_noob <- processData('Data/LFS_Mut/850k/schiffman',id_map_clin,"noob","beta","850")
saveRDS(data_schiffman_noob,'rds/Schiffman_Noob_beta.rds')

cat("[ Process 850k batch5]","\n")
data_batch5_raw <- processData('Data/LFS_Mut/850k/batch5',id_map_clin,"raw","beta","850")
data_batch5_noob <- processData('Data/LFS_Mut/850k/batch5',id_map_clin,"noob","beta","850")
saveRDS(data_batch5_noob,'rds/Batch5k_Noob_beta.rds')

cols <- Reduce(intersect, list(names(data_450k_noob),names(data_schiffman_noob),names(data_batch5_noob)))
data_450k_raw <- data_450k_raw[cols] ; data_schiffman_raw <- data_schiffman_raw[cols] ; data_batch5_raw <- data_batch5_raw[cols]
data_450k_noob <- data_450k_noob[cols] ; data_schiffman_noob <- data_schiffman_noob[cols] ; data_batch5_noob <- data_batch5_noob[cols]

cat("[ Combine data ]","\n")
dataRaw <- rbind(data_450k_raw,data_schiffman_raw,data_batch5_raw) ; saveRDS(dataRaw,'rds/Raw_beta.rds')
dataNoob <- rbind(data_450k_noob,data_schiffman_noob,data_batch5_noob) ; saveRDS(dataNoob,'rds/Noob_beta.rds')

ind <- 67
dataRaw <- readRDS('rds/Raw_beta.rds')
dataNoob <- readRDS('rds/Noob_beta.rds')

BRaw_cleaned <- dataRaw[(ind+1):length(dataRaw)]
BRaw_cleaned <- BRaw_cleaned[, colSums(is.na(BRaw_cleaned)) == 0]
B_cleaned <- dataNoob[(ind+1):length(dataNoob)]
B_cleaned <- B_cleaned[, colSums(is.na(B_cleaned)) == 0]
pheno <- dataNoob[1:ind]

########################## PLOT AND PROCESS RAW ##########################

cat("[ Plot PCs ]","\n")
indices <- which(is.na(pheno$cancer_diagnosis))

pc <- prcomp(as.matrix(BRaw_cleaned[-indices,]), scale = TRUE)
pc_clin <- cbind(pheno[-indices,],pc$x)
write.csv(pc_clin,paste0('Output/Raw_beta_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,"Raw_beta_PCA_summary.csv",outdir)
generate_pcplots(pc_clin,"Raw_beta",outdir)

########################## PLOT AND PROCESS NOOB ##########################

pc <- prcomp(as.matrix(B_cleaned[-indices,]), scale = TRUE)
pc_clin <- cbind(pheno[-indices,],pc$x)
write.csv(pc_clin,paste0('Output/Noob_beta_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,"Noob_beta_PCA_summary.csv",outdir)
generate_pcplots(pc_clin,"Noob_beta",outdir)

########################## REMOVE SEX CHROMOSOME PROBES ##########################

B_nosex <- remove_sex(B_cleaned, pheno, "Noob_beta")

pc <- prcomp(as.matrix(B_nosex[-indices,]), scale = TRUE)
pc_clin <- cbind(pheno[-indices,],pc$x)
write.csv(pc_clin,paste0('Output/Noob_beta_NoSex_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,"Noob_beta_NoSex_PCA_summary.csv",outdir)
generate_pcplots(pc_clin,"Noob_beta_NoSex",outdir)

cat("[ Plotting technical replicates ]","\n")

duplicated_data <- get_technicalreplicates(dataNoob)
plot_concordance(duplicated_data,"Noob_beta",outdir)

pc_beta <- data.frame(duplicated_data[(ind+1):length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:ind],pc$x)
write.csv(pc_clin,'Output/Noob_beta_TechnicalReplicates_PCA.csv',quote=F,row.names=F)

generate_pcsummary(pc_clin,'Noob_beta_TechnicalReplicates_PCA_summary.csv',outdir)
generate_pcplots(pc_clin,'Noob_beta_TechnicalReplicates',outdir)

