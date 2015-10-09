setwd('/home/yeung/projects/tissue-specificity')


# Run clustesr ------------------------------------------------------------

load("Robjs/bicmat.11_tiss_max_3.Robj")
n.centers <- 75
clusters <- kmeans(bic.mat, centers = n.centers)

save(clusters, file = "Robjs/bicmat.11_tiss_max_3.clusters.Robj")
