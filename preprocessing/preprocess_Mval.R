suppressMessages(require(bumphunter))
suppressMessages(require(minfi))
suppressMessages(require(dplyr))
suppressMessages(require(IlluminaHumanMethylation450kmanifest))
suppressMessages(require(IlluminaHumanMethylationEPICmanifest))
suppressMessages(require(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
suppressMessages(require(ggplot2))
suppressMessages(require(limma))
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

cat("[ Reading in methylation data]","\n")

cat("[ Process 450k ]","\n")
data_450k_raw <- processData('Data/LFS_Mut/450k',id_map_clin,"raw","M","450")
data_450k_noob <- processData('Data/LFS_Mut/450k',id_map_clin,"noob","M","450")
saveRDS(data_450k_noob,'rds/450k_Noob_Mval.rds')

cat("[ Process 850k Schiffman ]","\n")
data_schiffman_raw <- processData('Data/LFS_Mut/850k/schiffman',id_map_clin,"raw","M","850")
data_schiffman_noob <- processData('Data/LFS_Mut/850k/schiffman',id_map_clin,"noob","M","850")
saveRDS(data_schiffman_noob,'rds/Schiffman_Noob_Mval.rds')

cat("[ Process 850k batch5]","\n")
data_batch5_raw <- processData('Data/LFS_Mut/850k/batch5',id_map_clin,"raw","M","850")
data_batch5_noob <- processData('Data/LFS_Mut/850k/batch5',id_map_clin,"noob","M","850")
saveRDS(data_batch5_noob,'rds/Batch5k_Noob_Mval.rds')

cat("[ Combine data ]","\n")
cols <- Reduce(intersect, list(names(data_450k_noob),names(data_schiffman_noob),names(data_batch5_noob)))
data_450k_raw <- data_450k_raw[cols] ; data_schiffman_raw <- data_schiffman_raw[cols] ; data_batch5_raw <- data_batch5_raw[cols]
data_450k_noob <- data_450k_noob[cols] ; data_schiffman_noob <- data_schiffman_noob[cols] ; data_batch5_noob <- data_batch5_noob[cols]
dataRaw <- rbind(data_450k_raw,data_schiffman_raw,data_batch5_raw) ; saveRDS(dataRaw,'rds/Raw_Mval.rds')
dataNoob <- rbind(data_450k_noob,data_schiffman_noob,data_batch5_noob) ; saveRDS(dataNoob,'rds/Noob_Mval.rds')

#dataRaw <- readRDS('rds/Raw_Mval.rds') 
#dataNoob <- readRDS('rds/Noob_Mval.rds')

ind <- 67
MRaw_cleaned <- dataRaw[(ind+1):length(dataRaw)]
MRaw_cleaned <- MRaw_cleaned[, colSums(is.na(MRaw_cleaned)) == 0]
M_cleaned <- dataNoob[(ind+1):length(dataNoob)]
M_cleaned <- M_cleaned[, colSums(is.na(M_cleaned)) == 0]
pheno <- dataNoob[1:ind]


########################## PLOT AND PROCESS RAW ##########################

indices <- which(is.na(pheno$cancer_diagnosis))

pc <- prcomp(as.matrix(MRaw_cleaned[-indices,]),na.action=na.omit,scale=TRUE)
pc_clin <- cbind(pheno[-indices,], pc$x)
write.csv(pc_clin,paste0('Output/Raw_M_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,"Raw_M_PCA_summary.csv",outdir)
generate_pcplots(pc_clin,"Raw_M",outdir)

########################## PLOT AND PROCESS NOOB ##########################

pc <- prcomp(as.matrix(M_cleaned[-indices,]),scale = TRUE)
pc_clin <- cbind(pheno[-indices,],pc$x)
write.csv(pc_clin,paste0('Output/Noob_M_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,"Noob_M_PCA_summary.csv",outdir)
generate_pcplots(pc_clin,"Noob_M",outdir)

########################## REMOVE SEX CHROMOSOME PROBES ##########################

M_nosex <- remove_sex(M_cleaned, pheno, "Noob_M")
M_nosex <- M_nosex[,!apply(M_nosex, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]

pc <- prcomp(as.matrix(M_nosex[-indices,]), scale = TRUE)
pc_clin <- cbind(pheno[-indices,],pc$x)
write.csv(pc_clin,paste0('Output/Noob_M_NoSex_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,"Noob_M_NoSex_PCA_summary.csv",outdir)
generate_pcplots(pc_clin,"Noob_M_NoSex",outdir)

cat("[ Plotting technical replicates ]","\n")

duplicated_data <- get_technicalreplicates(dataM)
plot_concordance(duplicated_data,"Noob_M")

pc_beta <- data.frame(duplicated_data[(ind+1):length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc_beta <- Filter(function(x) sd(x) != 0, pc_beta)
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:int],pc$x)
write.csv(pc_clin,'Output/Noob_M_TechnicalReplicates_PCA.csv',quote=F,row.names=F)

generate_pcsummary(pc_clin,'Noob_M_TechnicalReplicates_PCA_summary.csv',outdir)
generate_pcplots(pc_clin,'Noob_M_TechnicalReplicates',outdir)

