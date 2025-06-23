library(maftools)
library(stringr)
source('custom_maftools.R')

dir <- "maftools_input/"

############################ Plot all pathogenic SNVs/indels in TP53 ############################
allmutsTable <- read.csv(paste0(dir,'lfs_meth_tp53_maf.csv'),stringsAsFactors = F)
p53_cnTable <- allmutsTable[c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification")]
p53_cnTable <- p53_cnTable[p53_cnTable$Variant_Classification == "Deletion",]
colnames(p53_cnTable) <- c("Gene","Sample_name","CN")
insdel_ids <- allmutsTable$Tumor_Sample_Barcode[allmutsTable$Variant_Classification=="Frame_Shift_Del_Ins"]
vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","Frame_Shift_Del_Ins",
                 "In_Frame_Ins", "Missense_Mutation","Deletion")
p53_maf <- read.maf(maf=allmutsTable,
                    cnTable=p53_cnTable,
                    clinicalData = paste0(dir,'lfs_meth_clinical.txt'),
                    vc_nonSyn = vc.nonSilent,
                    useAll = TRUE)
lollipopPlot_custom(maf = p53_maf, 
             gene = 'TP53', 
             showMutationRate = FALSE,
             refSeqID = 'NM_000546',
             AACol='aaChange',
             labelPos=c('282','175','273','220','245','248'))
plotmafSummary_custom(maf = p53_maf,  addStat = 'median',titvRaw = FALSE)

################################ Plot all pathogenic SNVs/indels ################################

#plot titv summary
laml.titv = titv(maf = p53_maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)

#oncoplot for top ten mutated genes.
oncoplot(maf = p53_maf, top = 10,clinicalFeatures=c('p53','cancer_diagnosis','family','gender'),sortByAnnotation = TRUE)

#oncogenic pathway analysis
OncogenicPathways(maf = p53_maf)
PlotOncogenicPathways(maf = p53_maf, pathways = "TP53")
