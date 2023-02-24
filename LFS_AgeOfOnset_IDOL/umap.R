suppressMessages(library(umap))

data <- readRDS("./valli_dataset.rds") # load
data_sm <- data[,45:ncol(data)]
data_sc_cen <- scale(data_sm, scale=TRUE, center=TRUE) # scale and center

one_umap <- function(df, n, dist){ #calculate umap and merge with clinical info and save
	u <- umap(df, n_neighbors=n, min_dist=dist) # umap
	df <- data.frame(x = u$layout[,1],
                   y = u$layout[,2])
	umap_feat <- merge(df, data[,1:44], by=0) # merge 
	filename = paste(umap, dist, sep="_")
	saveRDS(umap_feat, file = filename)
}

n_seq=c(0.2, 0.4, 0.6, 0.8)
for (x in n_seq){
	print(n_seq)
	one_umap(data_sc_cen, 15, x) # neighbours is 15, distance changes
}
