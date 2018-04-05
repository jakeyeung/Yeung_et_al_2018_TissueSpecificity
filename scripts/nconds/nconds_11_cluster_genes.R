# 2015-10-01
# Jake Yeung

setwd("/home/yeung/projects/tissue-specificity")

load("Robjs/bicmat.sub.Robj", verbose = T)
# d <- dist(t(bicmat.sub), method = "euclidean") # distance matrix
d <- dist(bicmat.sub, method = "euclidean") # distance matrix by gene
clustfit <- hclust(d, method="ward") 
# plot(clustfit) # display dendogram  # too many genes to plot
groups <- cutree(clustfit, k=40) # cut tree into 120 clusters
# draw dendogram with red borders around the 5 clusters 
# rect.hclust(clustfit, k=50, border="red")
# save(d, clustfit, file = "Robjs/nconds.cluster_2048_models.Robj")
save(d, clustfit, file = "Robjs/nconds.cluster_2048_models.by_gene.Robj")