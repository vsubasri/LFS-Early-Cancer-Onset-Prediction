suppressMessages(library(library(ggplot2))
suppressMessages(library(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))

# Working directory should be project root when called from shell scripts

parser <- ArgumentParser()
parser$add_argument("--datafile", action="store")
parser$add_argument("--id", action="store")
args <- parser$parse_args()

data <- readRDS(args$datafile)
beta <- data[45:length(data)]

## Plot PCA 
pc <- prcomp(as.matrix(beta), scale = TRUE,center=TRUE)
pc_clin <- cbind(pc$x,data[1:44])

## Find pcs correlated with age of sample collection 
coeff <- list(); pval <- list(); corr <- list()
for (i in 1:dim(pc$x)[1]) {
  lm_pc_age <- lm(scale(data$age_sample_collection) ~ scale(pc$x[, i]))
  pval[[i]] <- summary(lm_pc_age)$coefficients[2,4]
  coeff[[i]] <- summary(lm_pc_age)$coefficients[2,1]
  corr[[i]] <- cor.test(scale(data$age_sample_collection),scale(pc$x[, i]),method="spearman")$p.value
}

pca_age <- data.frame(lm_pval=do.call('rbind', pval), 
                      lm_coeff=do.call('rbind', coeff),
                      spear_pval=do.call('rbind', corr))
pca_age$rank <- rank(pca_age$spear_pval)
pca_age$pc <- paste0("PC",row.names(pca_age))

maxpc <- pca_age$pc[pca_age$rank == 1]
cat(paste0("[ Max age PC ] :",maxpc,"\n"))
maxpc2 <- pca_age$pc[pca_age$rank == 2]
cat(paste0("[ Second max age PC ] :",maxpc2,"\n"))

## Plot the pcs with the highest age of sample collection correlation

pdf(paste0("Plots/",args$id,".pdf"),width=9,height=7)
ggplot(pc_clin,aes_string(x=maxpc,y=maxpc2,color="age_sample_collection")) +
  geom_point() +
  theme_bw() + 
  scale_color_gradient(high="red",low="green")
suppressMessages(dev.off())

cat("Regress out PC most correlated with age",'\n')
##regress out pc most correlated to age
#https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com
XhatT <- pc$x[, !(colnames(pc$x) %in% c(maxpc))] %*% t(pc$rotation[, !(colnames(pc$rotation) %in% c(maxpc))])
beta_minusage <- t( XhatT * pc$scale + pc$center) %>% data.frame()
data_minusage <- cbind(data[1:44],t(beta_minusage))
saveRDS(data_minusage,paste0("rds/",args$id,"_MaxAgePCAdj",".rds"))

pc_minusage <- prcomp(as.matrix(t(beta_minusage)), scale = TRUE,center=TRUE)
pc_clin_minusage <- cbind(data[1:44],pc_minusage$x)

## Find pcs correlated with age of sample collection after correction
coeff_minusage <- list(); pval_minusage <- list(); corr_minusage <- list()
for (i in 1:dim(pc_minusage$x)[1]) {
  lm_pc_age <- lm(scale(data_minusage$age_sample_collection) ~ scale(pc_minusage$x[, i]))
  pval_minusage[[i]] <- summary(lm_pc_age)$coefficients[2,4]
  coeff_minusage [[i]] <- summary(lm_pc_age)$coefficients[2,1]
  corr_minusage[[i]] <- cor.test(scale(data_minusage$age_sample_collection),scale(pc_minusage$x[, i]),method="spearman")$p.value
}

pca_age_minusage <- data.frame(lm_pval_minusage=do.call('rbind', pval_minusage), 
                      lm_coeff_minuage=do.call('rbind', coeff_minusage),
                      spear_pval_minusage=do.call('rbind', corr_minusage))
pca_age_minusage$rank <- rank(pca_age_minusage$spear_pval)
pca_age_minusage$pc <- paste0("PC",row.names(pca_age_minusage))
write.csv(pca_age_minusage,paste0("Plots/",args$id,"_MaxAgePCAdj_PCA.csv"),quote=F,row.names=F)

maxpc_minusage <- pca_age_minusage$pc[pca_age_minusage$rank == 1]
maxpc2_minusage <- pca_age_minusage$pc[pca_age_minusage$rank == 2]

## Plot pcs with the highest age of sample collection correlation after correction

pdf(paste0("Plots/",args$id,"_MaxAgePCAdj",".pdf"),width=9,height=7)
ggplot(pc_clin_minusage,aes_string(x=maxpc_minusage,y=maxpc2_minusage,color="age_sample_collection")) +
  geom_point() +
  theme_bw() + 
  scale_color_gradient(high="red",low="green")
suppressMessages(dev.off())

