suppressMessages(library(minfi))
suppressMessages(library(FlowSorted.Blood.EPIC))
suppressMessages(library(IlluminaHumanMethylation450kmanifest))
suppressMessages(library(IlluminaHumanMethylationEPICmanifest))
suppressMessages(library(ExperimentHub))
#hub <- ExperimentHub()
#query(hub, "FlowSorted.Blood.EPIC")
#FlowSorted.Blood.EPIC <- hub[["EH1136"]]
#FlowSorted.Blood.EPIC

# find sample sheets
targets_450 <- read.metharray.sheet("data/Data/LFS_450", pattern = "SampleSheet.csv", ignore.case = TRUE, recursive = TRUE, verbose = TRUE)
targets_850 <- read.metharray.sheet("data/Data/LFS_850", pattern = "SampleSheet.csv", ignore.case = TRUE, recursive = TRUE, verbose = TRUE)

# remove rows where base name column is empty because samples listed in sample sheet but no idat files exist
targets_450<- targets_450[targets_450$Basename != "character(0)", ]
targets_850<- targets_850[targets_850$Basename != "character(0)", ]
# removed 132 samples (most from cases_toronto)
saveRDS(targets_450, file = "450.rds")
saveRDS(targets_850, file = "850.rds")

#targets_450 <- readRDS("450.rds") # uncomment if ran and saved above
#targets_850 <- readRDS("850.rds")

# create RGChannelSet object for IDOL
RGsetTargets_450 <- read.metharray.exp(targets = targets_450) 
RGsetTargets_850 <- read.metharray.exp(targets = targets_850, force=TRUE) 

# combine RGset arrays
# out type can also be "IlluminaHumanMethylationEPIC"
# converting to 450k tested 
# removing 7% of probes in 450k dataset did not change cell type proportions
# correlation of cell type proportions before and after removing 7% of probes was >99% 
# 850k contains 93% pf 450k probes
RGsetTargets <- combineArrays(RGsetTargets_450, RGsetTargets_850, outType = "IlluminaHumanMethylation450k", verbose = TRUE)
RGsetTargets_EPIC <- combineArrays(RGsetTargets_450, RGsetTargets_850, outType = "IlluminaHumanMethylationEPIC", verbose = TRUE)
saveRDS(RGsetTargets_EPIC, file = "targets_EPIC.rds")

# load epic targets
#RGsetTargets_EPIC<- readRDS("targets_EPIC.rds") # uncomment if already ran and saved above

# estimate immune cell abundance
# noob preprocess recommended by IDOL
# can change IDOL to “both” and use cell types "Bcell", "CD4T", "CD8T", "Eos", "Gran", "Mono", "Neu", "NK" when change to “both”
data(IDOLOptimizedCpGs450klegacy)
data(IDOLOptimizedCpGs)
#count_both<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood", processMethod = "preprocessNoob", probeSelect = "IDOL", cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), referencePlatform="IlluminaHumanMethylation450k", referenceset=NULL, IDOLOptimizedCpGs=IDOLOptimizedCpGs450klegacy) # if using RGset with 450k array
count_both<-estimateCellCounts2(RGsetTargets_EPIC, compositeCellType = "Blood", processMethod = "preprocessNoob", probeSelect = "IDOL", cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), referencePlatform="IlluminaHumanMethylationEPIC", referenceset=NULL, IDOLOptimizedCpGs=IDOLOptimizedCpGs) # if using RGset with EPIC array

saveRDS(count_both, file = "count_EPIC.rds")

## extended cell types
#prop_ext<-estimateCellCounts2(RGsetTargets_EPIC, compositeCellType = "BloodExtended", processMethod = "preprocessNoob", probeSelect = "IDOL", cellTypes = c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg"), referencePlatform="IlluminaHumanMethylationEPIC", referenceset=NULL, IDOLOptimizedCpGs=IDOLOptimizedCpGs)
##CustomCpGs =if(RGsetTargets@annotation[1]=="IlluminaHumanMethylationEPIC"){IDOLOptimizedCpGsBloodExtended else{IDOLOptimizedCpGsBloodExtended450k})
##perc_ext<-round(prop_ext$prop*100,1)
#
#saveRDS(prop_ext, file = "count_EPIC_ext.rds")
