
metadata <- read.csv('data/GSE42861_metadata.csv', header=F)
filtcols <- sub("\\:.*", "",metadata$V2)
metadata <- t(metadata)
colnames(metadata) <- metadata[1,]
metadata <- as.data.frame(metadata[-1,])
start <- ## starting index
end <- ## ending index
colnames(metadata)[start:end] <- filtcols[start:end] 

metadata[] <- lapply(metadata, function(x) sub(".*:", "", as.character(x)))
write.table(metadata,'data/GSE42861_metadata.txt', quote=F, sep='\t')
