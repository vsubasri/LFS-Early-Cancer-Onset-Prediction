library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(stringr)

features <- read.csv('data/features.txt',sep='\t')
ann450k <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
ann450k$region <- do.call(rbind,str_split(ann450k$UCSC_RefGene_Group,';'))[,1]
ann450k_3utr <- ann450k[ann450k$region=="3'UTR",]
ann450k_3utr <- data.frame(ann450k_3utr %>% group_by(chr) %>% summarize(n=n()))

probe_loc <- data.frame(probe = features$probe,
                chr =ann450k$chr[match(features$probe,ann450k$Name)])
probe_loc$chr <- factor(probe_loc$chr,levels=paste0("chr",c((1:22),"X")))
probe_summ <- probe_loc %>% group_by(chr) %>% summarize(n=n())
probe_summ$n_3utr <- ann450k_3utr$n[match(probe_summ$chr,probe_summ$chr)]
probe_summ$enrichment <- probe_summ$n/probe_summ$n_3utr

ggplot(probe_summ,aes(x=chr,y=enrichment)) + 
  geom_bar(color="black",stat="identity") + 
  theme_classic() +
  theme(text = element_text(size=16,family="Helvetica Neue"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Chromosome",y="# of Model Probes /\n# of Probes in 3'UTR")

