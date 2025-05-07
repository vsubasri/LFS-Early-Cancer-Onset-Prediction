require(e1071)
require(xgboost)
require(shapviz)
require(stats)
require(umap)
require(tidyverse)
require(dplyr)
library(kernelshap)
library(patchwork) 

model_name <- "svmRadial"
model_path <- paste0('checkpoint/NoobCorrected_beta_ProjPC2Adj_lfs_3UTR_tp53_72_svmRadial_ageofonset_results.rds')
data_path <- paste0("data/data_train.rds")
feature_path <- paste0("data/features.txt")

# load model
all <- readRDS(model_path)
train_x_y <- all[[1]]$trainingData # training set with x and y
names(train_x_y)[names(train_x_y) == '.outcome'] <- 'cancer_before_six' # outcome column in training set
training_data <- subset(all[[1]]$trainingData, select=-c(.outcome)) # training dataset drop outcome 
train <- xgb.DMatrix(as.matrix(training_data)) # training set format for xgboost shap
if (model_name == "xgboost") {
  xgb.save(all[[1]]$finalModel, 'ageofonset.model') # save model in diff format
  model <- xgb.load('ageofonset.model') # load model
} else {
  model <- all[[1]]
}
test_data_pred <- all[[2]] # test set predictions

# load data
data <- readRDS(data_path)
# load model probes and subset data
probes <- read.csv(feature_path, sep="\t", row.names=NULL)
colnames(probes) <- c("probe", "gene")
probes_model <- probes$probe # probe list to filter dataset
#data_probes <- select(data, probes_model) # filter data by model probes
#data_filt <- cbind(data[1:ind], data_probes) # data with clinical and probes in model

### SHAP ###
if (model_name %in% c("xgboost", "gbm")) {
  shap_xg <- predict(model, train, predcontrib = TRUE) # xgboost shap values
  shap_model <- subset(shap_xg, select=-c(BIAS)) # remove extra column
} else if (model_name %in% c("svmRadial","svmLinear","enet")) {
  pred_fun <- function(mod, X) predict(model,X,type="prob")[2]
  shap_model <- kernelshap(model, training_data, bg_X = training_data, pred_fun = pred_fun) 
} else {
  print("Model not supported")
}

#saveRDS(shap_model,'shap_model.rds')
shap_model <- readRDS(paste0(dir,'data/shap_model.rds'))
validation <- readRDS(paste0(dir,"data/data_val.rds"))
validation_data <- validation[colnames(training_data)]
test <- readRDS(paste0(dir,"data/data_test.rds"))
positive_rows <- which(test[["ageofonset"]] < 72)
test_data <- test[colnames(training_data)]
shp_val <- shapviz(shap_model, X = training_data, X_pred = validation_data)
shp_test <- shapviz(shap_model, X = training_data, X_pred = test_data)

# visualize dataset
# left side is cancer after 6
# right side is cancer before 6
sv_importance(shp_val, kind = "both", alpha=0.5, max_display = 40, show_numbers = TRUE) + theme_minimal()
sv_importance(shp_val, kind = "bar", max_display=130) + theme_minimal() 

sv_importance(shp_test, kind = "both", alpha=0.5, max_display = 40, show_numbers = TRUE) + theme_minimal()
sv_importance(shp_test, kind = "both", size=2, max_display=1000) + theme_minimal()+ scale_x_continuous(position = "top")

feats <- data.frame(shap=apply(shp_test$S, 2, function(x) mean(abs(x)))) 
feats$feature <- rownames(feats)
ggplot(feats[feats$feature %in% c("dbd","tad","od","nterm","missense","splice","nonsense","frameshift","deletion"),],aes(reorder(feature, shap),shap)) +
  geom_bar(stat = "identity") + 
  labs(x= "TP53 Mutation Feature", y = "Feature Importance") +
  theme_classic() +
  theme(text = element_text(size=16,family="Helvetica Neue"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() 

#t <- table(train_x_y$cancer_before_six, train_x_y$cancer_atdraw)
#rownames(t) <- c("cancer after 6", "cancer before 6") # table shows all patients without cancer at draw have cancer after 6 so good predictor

# visualize feature
sv_dependence(shp_test, v="LEF1")
sv_dependence(shp_test, v="BTBD7", "auto")
sv_dependence(shp_test, v="MSL3", "auto")
sv_dependence(shp_test, v="TSR1", "auto")
sv_dependence(shp_test, v="RET", "auto")

# visualize one datapoint
sv_waterfall(shp_test,row_id =positive_rows[8]) +
  theme(axis.text = element_text(size = 11))
a <- sv_force(shp_test, row_id = positive_rows[1])
b <- sv_force(shp_test, row_id = positive_rows[2])
c <- sv_force(shp_test, row_id = positive_rows[3])
d <- sv_force(shp_test, row_id = positive_rows[4])
e <- sv_force(shp_test, row_id = positive_rows[5])
f <- sv_force(shp_test, row_id = positive_rows[6])
g <- sv_force(shp_test, row_id = positive_rows[7])
h <- sv_force(shp_test, row_id = positive_rows[8]) ## incorrect
i <- sv_force(shp_test, row_id = positive_rows[9])
j <- sv_force(shp_test, row_id = positive_rows[10])
ggarrange(a,b,c,d,e,f,g,h,i,j, nrow=5,ncol=2,labels = na.omit(test$ids[test$ageofonset< 72]))

# IL15RA and outcome
ggplot(train_x_y, aes(x=TSR1, y=cancer_before_six)) +
  geom_point(alpha=0.8)

### PCA ###
# with all probes on all samples
ind <- 68
pca <- prcomp(data[,ind:ncol(data)], scale=TRUE, center=TRUE) # exclude systemic treatment because all 0 and outcome column
summary(pca) # need 53 components to describe 90% of the variance
# 2 PC only rep 37% variance

# coloured by age of onset
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$ageofonset)) +
  geom_point(alpha=0.8, size=2) +
  labs(color='Age of onset') +
  theme_classic(base_size=18)

# coloured by cancer_diagnosis 
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$cancer_diagnosis)) +
  geom_point(alpha=0.8, size=2) +
  labs(color='Cancer diagnosis') +
  theme_classic(base_size=18)

# coloured by tissue type 
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$tissue_type)) +
  geom_point(alpha=0.8, size=2) +
  labs(color='Cancer diagnosis') +
  theme_classic(base_size=18)

# coloured by cancer_atdraw
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$cancer_atdraw)) +
  geom_point(alpha=0.8, size=2) +
  labs(color='Cancer at draw') +
  theme_classic(base_size=18)

# coloured by project
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$Project)) +
  geom_point(alpha=0.8, size=2) +
  labs(color='Project') +
  theme_classic(base_size=18)

### UMAP ###
data_sm <- data[,(ind+1):ncol(data)]
data_sc_cen <- scale(data_sm, scale=TRUE, center=TRUE)
u <- umap(data_sc_cen) # tried diff params
df <- data.frame(x = u$layout[,1],
                 y = u$layout[,2])
umap_feat <- merge(df, data[,1:ind], by=0)

# coloured by age of onset
ggplot(umap_feat, aes(x, y, colour = ageofonset)) +
  geom_point(size=2) +
  theme_classic()

# coloured by cancer_diagnosis 
ggplot(umap_feat, aes(x, y, colour = cancer_diagnosis)) +
  geom_point(size=2) +
  theme_classic()

# coloured by tissue type
ggplot(umap_feat, aes(x, y, colour = tissue_type)) +
  geom_point(size=2) +
  theme_classic()

# coloured by cancer_atdraw
ggplot(umap_feat, aes(x, y, colour = cancer_atdraw)) +
  geom_point(size=2) +
  theme_classic()

# coloured by project
ggplot(umap_feat, aes(x, y, colour = Project)) +
  geom_point(size=2) +
  theme_classic()


