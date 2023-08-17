library(xgboost)
library(shapviz)
library(stats)
library(umap)
library(tidyverse)
library(dplyr)

model_path <- ##INSERT HERE
data_path <- ##INSERT HERE
feature_path <- ##INSERT HERE
ind <- 67

# load model
all <- readRDS(model_path)
train_x_y <- all[[1]]$trainingData # training set with x and y
names(train_x_y)[names(train_x_y) == '.outcome'] <- 'cancer_before_six' # outcome column in training set
training_data <- subset(all[[1]]$trainingData, select=-c(.outcome)) # training dataset drop outcome 
train <- xgb.DMatrix(as.matrix(training_data)) # training set format for xgboost shap
xgb.save(all[[1]]$finalModel, 'ageofonset.model') # save model in diff format
model <- xgb.load('ageofonset.model') # load model
test_data_pred <- all[[2]] # test set predictions

# load data
data <- readRDS(data_path)
# load model probes and subset data
probes <- read.csv(feature_path, sep="\t", row.names=NULL) 
colnames(probes) <- c("probe", "gene")
probes_model <- probes$probe # probe list to filter dataset
data_probes <- select(data, probes_model) # filter data by model probes
data_filt <- cbind(data[1:ind], data_probes) # data with clinical and probes in model

### SHAP ###
levels(train_x_y$cancer_before_six) # levels for shap
all[[1]]$finalModel$obsLevels # levels
shap_xg <- predict(model, train, predcontrib = TRUE) # xgboost shap values
shap_xgb <- subset(shap_xg, select=-c(BIAS)) # remove extra column

# shapviz
shp <- shapviz(model, X_pred = data.matrix(training_data), X = training_data) # calculate with shapviz
shp <- shapviz(shap_xgb, X = training_data) # use shap values from xgboost

# visualize dataset
# left side is cancer after 6
# right side is cancer before 6
sv_importance(shp, kind = "beeswarm", alpha=0.5, show_other = FALSE)
sv_importance(shp, kind = "bar", show_other = FALSE)
t <- table(train_x_y$cancer_before_six, train_x_y$cancer_atdraw)
rownames(t) <- c("cancer after 6", "cancer before 6") # table shows all patients without cancer at draw have cancer after 6 so good predictor

# visualize feature
sv_dependence(shp, v="cancer_atdraw")
sv_dependence(shp, v="IL15RA", "auto")
sv_dependence(shp, v="SEMA7A", "auto")

# visualize one datapoint
sv_waterfall(shp, row_id = 1)
sv_force(shp, row_id = 1)

# IL15RA and outcome
ggplot(train_x_y, aes(x=IL15RA, y=cancer_before_six)) + 
  geom_point(alpha=0.8)



### PCA ###
# with all probes on all samples
pca <- prcomp(data_filt[,45:ncol(data_filt)], scale=TRUE, center=TRUE) # exclude systemic treatment because all 0 and outcome column
summary(pca) # need 53 components to describe 90% of the variance
# 2 PC only rep 37% variance

# coloured by age of onset
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$ageofonset)) + 
  geom_point(alpha=0.8, size=5) + 
  labs(color='Age of onset') +
  theme_classic(base_size=18) 

# coloured by cancer_diagnosis 
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$cancer_diagnosis)) + 
  geom_point(alpha=0.8, size=5) + 
  labs(color='Cancer diagnosis') +
  theme_classic(base_size=18) 

# coloured by cancer_atdraw
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$cancer_atdraw)) + 
  geom_point(alpha=0.8, size=5) + 
  labs(color='Cancer at draw') +
  theme_classic(base_size=18) 

# coloured by project
ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2, color=data$Project)) + 
  geom_point(alpha=0.8, size=5) + 
  labs(color='Project') +
  theme_classic(base_size=18) 



### UMAP ###
data_sm <- data_filt[,(ind+1):ncol(data_filt)]
data_sc_cen <- scale(data_sm, scale=TRUE, center=TRUE)
u <- umap(data_sc_cen, n_neighbours=15, min_dist=0.6) # tried diff params
df <- data.frame(x = u$layout[,1],
                 y = u$layout[,2])
umap_feat <- merge(df, data_filt[,1:ind], by=0)

# coloured by age of onset
ggplot(umap_feat, aes(x, y, colour = ageofonset)) +
  geom_point(size=4) +
  theme_classic() 

# coloured by cancer_diagnosis 
ggplot(umap_feat, aes(x, y, colour = cancer_diagnosis)) +
  geom_point(size=4) +
  theme_classic()

# coloured by cancer_atdraw
ggplot(umap_feat, aes(x, y, colour = cancer_atdraw)) +
  geom_point(size=4) +
  theme_classic()

# coloured by project
ggplot(umap_feat, aes(x, y, colour = Project)) +
  geom_point(size=4) +
  theme_classic()
