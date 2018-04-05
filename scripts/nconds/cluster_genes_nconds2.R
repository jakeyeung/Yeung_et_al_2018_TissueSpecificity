setwd('/home/yeung/projects/tissue-specificity')

library(reshape2)
library(dplyr)
# Run clustesr ------------------------------------------------------------

n.centers <- 65
# inpath <- "Robjs/bicmat.11_tiss_max_3.clusters.top5.bug_fixed.Robj"
outpath <- paste0("Robjs/bicmat.11_tiss_max_3.clusters.top5.bug_fixed.clusters.", n.centers, ".Robj")
load("/home/yeung/projects/nconds_results/fits_long.11_tiss_3_max.weight_raw.Robj", verbose=T)  # fits.long
fits.long <- subset(fits.long, weight > 0)

fits.long <- fits.long %>%
  group_by(gene) %>%
  arrange(desc(weight)) %>%
  do(head(., n = 10))

print(head(data.frame(fits.long, n = 100)))

mats.all <- dcast(fits.long, gene ~ model, value.var = "weight", fill = 0)
save(mats.all, file = "Robjs/mats.all.cluster.top10.Robj")

# fix gene names
rownames(mats.all) <- mats.all$gene; mats.all$gene <- NULL

print(mats.all[1:5, 1:5])
mats.all <- mats.all[complete.cases(mats.all), ]
clusters <- kmeans(t(mats.all), centers = n.centers)

save(clusters, file = outpath)
