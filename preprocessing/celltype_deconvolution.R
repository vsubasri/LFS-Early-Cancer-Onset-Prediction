library(minfi)
suppressMessages(library(gplots))
suppressMessages(library(mlogit))
suppressMessages(library(dplyr))

# Working directory should be project root when called from shell scripts

source('preprocessing/tumor_type_functions.R')

################################ READ IN INPUT ID FILES #################################

cat("[ Reading in input sample sheets ]","\n")

ind <- 
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
# Estimate Cell Counts
###################################################################################################################

cat("[ Estimating cell counts ]","\n")

if(!file.exists('rds/celltype_props.rds')){

        rgCasesT <- read.metharray.exp('Data/LFS_Mut/450k/cases_toronto', recursive = T)
        celltype_propsT <- estimateCellCounts(rgCasesT)

        rgCasesM <- read.metharray.exp('Data/LFS_Mut/450k/cases_montreal', recursive = T)
        celltype_propsM <- estimateCellCounts(rgCasesM)

        rgCasesUtah <- read.metharray.exp('Data/LFS_Mut/850k/schiffman', recursive = T)
        celltype_propsUtah <- estimateCellCounts(rgCasesUtah)

        rgCasesbatch5 <- read.metharray.exp('Data/LFS_Mut/850k/batch5', recursive = T)
        celltype_propsbatch5 <- estimateCellCounts(rgCasesbatch5)

        celltype_props <- rbind(celltype_propsM, celltype_propsT,celltype_propsbatch5,celltype_propsUtah)
        saveRDS(celltype_props, 'rds/celltype_props.rds')

} else {

        celltype_props <- readRDS('rds/celltype_props.rds')
}


###################################################################################################################
# Calculate NLR
###################################################################################################################

granulocytes = celltype_props[, "Gran"]
lymphocytes = celltype_props[, c("Bcell","CD8T","CD4T","NK")]

nlr = data.frame(nlr = granulocytes/rowSums(lymphocytes))
nlr <- transform(merge(nlr,data.frame(celltype_props),by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

###################################################################################################################
# Regress out NLR
###################################################################################################################

clin <- read.csv('clinical_two.csv', stringsAsFactors = F, fileEncoding="UTF-8-BOM")
nlr_T_M <- findIds(nlr[!rownames(nlr) %in% c(rownames(id_map_val),rownames(id_map_con)),], id_map) %>% getIdName() %>% cleanIds()
nlr_Val <- transform(merge(nlr[rownames(nlr) %in% rownames(id_map_val),],id_map_val[c("Sample.ID","Sentrix.Barcode")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
nlr_Con <- transform(merge(nlr[rownames(nlr) %in% rownames(id_map_con),],id_map_con[c("Sample_Name","Sentrix_ID")],by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

colnames(nlr_Val)[8:9] <- colnames(nlr_T_M)[8:9]
colnames(nlr_Con)[8:9] <- colnames(nlr_T_M)[8:9]

nlr <- rbind(nlr_T_M,nlr_Val,nlr_Con)

nlr$methid <- rownames(nlr)
nlr_clin<- merge(nlr,clin[!duplicated(clin$blood_dna_malkin_lab_),], by.x="ids", by.y="blood_dna_malkin_lab_",all.x=T) 
nlr_clin$cancerstatus <- ifelse(nlr_clin$cancer_diagnosis_diagnoses == "Unaffected", 0, 1)

celltypes= c("CD8T","CD4T","NK","Bcell","Mono","Gran")

nlr_log<- nlr_clin[!is.na(nlr_clin$cancerstatus),]
nlr_log[celltypes] <- apply(nlr_log[celltypes],2,scale)

####### Evaluating immune contexture under mutant status only #######
nlr_tumortype <- nlr_log_mut <- nlr_log[nlr_log$p53_germline == "Mut" & !is.na(nlr_log$p53_germline),]
fmlt_int <- as.formula(paste0("cancerstatus ~ CD8T + CD4T + NK + Bcell + Mono + Gran"))
fmlt_noint <- as.formula(paste0("cancerstatus ~ 0 + CD8T + CD4T + NK + Bcell + Mono + Gran"))
cancertype.log <- glm(fmlt_int, family=binomial(link='logit'), data=nlr_log,trace=FALSE)
celltype_assoc <- data.frame(coef(summary(cancertype.log)))
celltype_assoc$celltype <- rownames(celltype_assoc)
celltype_assoc <- celltype_assoc[celltype_assoc$celltype != "(Intercept)",]


ggplot(celltype_assoc,aes(x=celltype,y=Estimate)) + 
  geom_bar(stat="identity",fill="red",alpha=0.6) + 
  theme_bw() + 
  labs(x ="Cell type", y = "Coefficient")
write.csv(celltype_assoc, "~/research/lfs_methylation/Output/celltype_assoc.csv", quote=F)

###################################################################################################################
# Using mlogit
###################################################################################################################

nlr_tumortype$cancer_diagnosis_diagnoses <- ifelse(grepl("DCIS|IDC|hylloides|mammary", nlr_tumortype$cancer_diagnosis_diagnoses), "Breast ca", nlr_tumortype$cancer_diagnosis_diagnoses)
nlr_tumortype$cancer_diagnosis_diagnoses <- ifelse(grepl("ARMS|ERMS|RMS", nlr_tumortype$cancer_diagnosis_diagnoses), "RMS", nlr_tumortype$cancer_diagnosis_diagnoses)
nlr_tumortype <- nlr_tumortype[nlr_tumortype$cancer_diagnosis_diagnoses %in% c("ACC","Breast ca","CPC","OS","RMS","Unaffected"),]
nlr_mlogit <- mlogit.data(nlr_tumortype, choice="cancer_diagnosis_diagnoses",shape="wide",id.var="id",varying = NULL)
mform_int <- mFormula(formula(paste("cancer_diagnosis_diagnoses ~ 0 | CD8T + CD4T + NK + Bcell + Mono + Gran")))
mform_noint <- mFormula(formula(paste("cancer_diagnosis_diagnoses ~ 0 | 0 + CD8T + CD4T + NK + Bcell + Mono + Gran")))
mlogit.log <- mlogit(mform_int, reflevel="Unaffected", data=nlr_mlogit)
mlogit_table <- data.frame(summary(mlogit.log)$CoefTable)
mlogit_table[c("Tumortype","Celltype")] <- do.call(rbind,str_split(rownames(mlogit_table),":"))
write.csv(mlogit_table, "~/research/lfs_methylation/Output/mlogit_table_tumortype_mut_int_removedMultipleSamePatient.csv", quote=F)

###################################################################################################################
# Clustering and Heatmap
###################################################################################################################

library(ComplexHeatmap)

celltype_plot <- melt(nlr_tumortype,id.vars="cancer_diagnosis_diagnoses",measure.vars = celltypes)
celltype_plot$cancerstatus  <- ifelse(celltype_plot$cancer_diagnosis_diagnoses == "Unaffected","Unaffected","Cancer")
mlogit_plot <- mlogit_table[mlogit_table$Celltype != "(intercept)",]

label.df <- mlogit_plot[c("Celltype","Tumortype","Pr...z..")]
label.df$value <- 7
label.df$label <- ifelse(label.df$Pr...z..<0.05, "*","")
colnames(label.df) <- c("cancer_diagnosis_diagnoses","variable","pvalue","value","label")

ggplot(celltype_plot,aes(x=variable,y=value,fill=cancer_diagnosis_diagnoses)) + 
  geom_boxplot() + 
  theme_bw() + 
  labs(x ="Cell type", y = "Normalized cell proportion", fill="Tumor type") 
  + geom_text(data = label.df, aes(x=variable,y=value,group=cancer_diagnosis_diagnoses,label = label),position=position_dodge(width=0.9))


ggplot(mlogit_plot,aes(x=Celltype,y=Estimate,fill=Tumortype)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  theme_bw() +
  geom_errorbar(aes(ymin=Estimate-Std..Error, ymax=Estimate+Std..Error), width=.2,
                position=position_dodge(.9)) 
  

nlr_clin_mut <- nlr_clin[grepl("Mut",nlr_clin$p53_germline),]
celltype_props_mut <- celltype_props[rownames(celltype_props) %in% nlr_clin_mut$methid,]
rownames(celltype_props_mut) <- make.unique(paste0(nlr_clin$ids[match(rownames(celltype_props_mut),nlr_clin$methid)], "(",
                                  nlr_clin$cancer_diagnosis_diagnoses[match(rownames(celltype_props_mut),nlr_clin$methid)],")"))
labels <- rownames(celltype_props_mut)
celltype_props_mut <- apply(celltype_props_mut, 2, scale)
rownames(celltype_props_mut) <- labels

dend <- dist(celltype_props_mut, method = 'manhattan') %>%  
  hclust(method = "complete") %>% 
  as.dendrogram %>%
  hang.dendrogram(hang_height=0.1) %>%
  set("labels_cex", 0.8) %>%
  color_branches(k=6) 

dend <- rotate(dend, 1:nleaves(dend))

par(mar = c(3,3,2,7))

#plot clustering 

plot(dend, 
     main = "Clustered Methylation of Cancer Types", 
     horiz =  TRUE,  nodePar = list(cex = .007)) 

# plot heatmap 
par(mar=c(5.1, 4.1, 0.5, 0.5))
hmcols <- colorpanel(2750,"blue","black", "yellow")
heatmap.2(as.matrix(celltype_props_mut),  
          srtCol = 60,
          key = TRUE,
          lhei=c(1,4), lwid=c(1,5), keysize=0.5, key.par = list(cex=0.25),
          trace="none",
          dendrogram = "row",
          Rowv = dend,
          Colv = "NA", 
          col=hmcols,
          cexCol=1.0,
          cexRow=0.8,
          scale='row',
          density.info="histogram"
)    



