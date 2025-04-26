library(methylclock)
library(ggplot2)
library(reshape2)
library(corrplot)
library(ggpubr)

dir <- "/Users/vallijahsubasri/Documents/lfs_ageofonset/"
lfs <- readRDS(paste0(dir,'data/Noob_beta.rds'))
lfs$agesamplecollection <- lfs$agesamplecollection/12
lfs$ageofonset <- lfs$ageofonset/12
lfs <- lfs[!is.na(lfs$agesamplecollection),]
meth <- t(lfs[68:length(lfs)])
colnames(meth) <- lfs$SentrixID
cpgs.missing <- checkClocks(meth)
methage <- DNAmAge(meth, clocks=c("Horvath", "Hannum", "Levine","PedBE", "Wu", "TL", "BLUP","EN"),age=lfs$agesamplecollection)
methage$ageofonset <- lfs$ageofonset
methage$earlyonset <- ifelse(is.na(methage$ageofonset),"Cancer-Free",ifelse(methage$ageofonset < 6,"Early Onset","Late Onset"))

ggplot(methage,aes(age,ageAcc.TL,color=ageofonset)) +
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  theme_minimal()
canceronly <- methage[!is.na(methage$ageofonset),]
res <- round(cor(canceronly[3:length(canceronly)-1]),2)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

methage_plot <- melt(methage,id=c("id","age","ageofonset","earlyonset"))
methage_plot_acc <- methage_plot[grepl("Acc",methage_plot$variable),]
ggplot(methage_plot_acc,aes(ageofonset,value,color=ageofonset)) +
  geom_point() + 
  facet_wrap(~variable,ncol=6,scales = "free") + 
  geom_point(data=methage_plot_acc[is.na(methage_plot_acc$ageofonset),],aes(age,value),color="grey",size=1) +
  geom_smooth(method='lm',color="red",width=0.1,size=1) + 
  stat_cor(size=2)+
  theme_minimal() 

ggplot(methage_plot_acc,aes(age,value,color=ageofonset)) +
  geom_point() + 
  facet_wrap(~variable,ncol=6,scales = "free") + 
  geom_smooth(method='lm',color="red",width=0.1,size=1) + 
  stat_cor(size=2)+
  theme_minimal() 

methage_plot_u6 <- methage_plot[methage_plot$age < 6,]
ggplot(methage_plot_u6,aes(earlyonset,value)) +
  geom_boxplot() + 
  geom_point(aes(color=age))+
  facet_wrap(~variable,ncol=8,scales = "free") + 
  stat_compare_means(method.args = list(alternative = "less")) +
  theme_minimal() 

ggplot(methage_plot,aes(earlyonset,value)) +
  geom_boxplot() + 
  geom_point(aes(color=age))+
  facet_wrap(~variable,ncol=8,scales = "free") + 
  stat_compare_means(method.args = list(alternative = "less")) +
  theme_minimal() 
