setwd('/home/yeung/projects/tissue-specificity')


# Run clustesr ------------------------------------------------------------

n.centers <- 500
inpath <- "Robjs/bicmat.11_tiss_max_3.clusters.top5.bug_fixed.Robj"
outpath <- paste0("Robjs/bicmat.11_tiss_max_3.clusters.top5.bug_fixed.clusters.", n.centers, ".Robj")
load(inpath)

# fix gene names
rownames(mats.all) <- mats.all$gene; mats.all$gene <- NULL

print(mats.all[1:5, 1:5])
mats.all <- mats.all[complete.cases(mats.all), ]
clusters <- kmeans(t(mats.all), centers = n.centers)

save(clusters, file = outpath)
