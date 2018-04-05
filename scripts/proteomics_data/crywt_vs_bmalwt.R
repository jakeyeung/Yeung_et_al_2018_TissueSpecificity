# 2016-11-24
# crywt_vs_bmalwt.R
# Show correlation between crywt and bmalwt in proteomics

rm(list=ls())

library(ggplot2)
library(reshape2)

source("scripts/functions/ProteomicsFunctions.R")

prot <- LoadProteomicsData(as.long=FALSE)
prot.long <- LoadProteomicsData()


# Compare WT vs BmalWT ----------------------------------------------------

jgenos <- c("WT", "BmalWT")
jtimes <- unique(subset(prot.long, geno == "BmalWT")$time)
jtimes <- jtimes[which(jtimes == 18)]
prot.sub <- subset(prot.long, geno %in% jgenos & time %in% jtimes & gene != "")
prot.crybmal <- dcast(prot.sub, formula = gene + Protein.IDs + time ~ geno, value.var = "rel.abund")
prot.crybmal <- prot.crybmal[complete.cases(prot.crybmal), ]
prot.crybmal$time <- as.factor(prot.crybmal$time)
m <- ggplot(prot.crybmal, aes_string(x = jgenos[[2]], y = jgenos[[1]])) + geom_point(alpha = 0.2) + theme_bw() + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + facet_wrap(~time)
jcor <- cor(prot.crybmal[[jgenos[[1]]]], prot.crybmal[[jgenos[[2]]]])
print(m + ggtitle(paste("Correlation", jcor)))

# ggpairs(subset(prot.crybmal, select = c(-gene, -Protein.IDs)), colour='time', alpha=0.4)


# Compare replicates ------------------------------------------------------

jgenos <- c("WT")
jtimes <- unique(subset(prot.long, geno == "WT")$time)
prot.sub <- subset(prot.long, geno %in% jgenos & time %in% jtimes & gene != "")
prot.sub$time <- paste("ZT", prot.sub$time, sep = "")
prot.cry <- dcast(prot.sub, formula = gene + Protein.IDs ~ time, value.var = "rel.abund")
prot.cry <- prot.cry[complete.cases(prot.cry), ]
m <- ggplot(prot.cry, aes(x = ZT0, y = ZT24)) + geom_point(alpha = 0.2) + theme_bw() + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)


# Compare two timepoints in BmalWT ----------------------------------------

jgenos <- c("BmalWT")
jtimes <- unique(subset(prot.long, geno == jgenos)$time)
prot.sub <- subset(prot.long, geno %in% jgenos & time %in% jtimes & gene != "")
prot.sub$time <- paste("ZT", prot.sub$time, sep = "")
prot.bmalwt <- dcast(prot.sub, formula = gene + Protein.IDs ~ time, value.var = "rel.abund")
prot.bmalwt <- prot.bmalwt[complete.cases(prot.bmalwt), ]
m <- ggplot(prot.bmalwt, aes(x = ZT0, y = ZT6)) + geom_point(alpha = 0.2) + theme_bw() + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ggtitle("ZT0 vs ZT6 BmalWT")
print(m)


# PCA of four samples -----------------------------------------------------

jgenos <- c("BmalWT")
jtimes <- unique(subset(prot.long, geno == jgenos)$time)
prot.sub <- subset(prot.long, geno %in% jgenos & time %in% jtimes & gene != "")
prot.sub$time <- paste("ZT", prot.sub$time, sep = "")
prot.bmalwt <- dcast(prot.sub, formula = gene + Protein.IDs ~ time, value.var = "rel.abund")
prot.bmalwt <- prot.bmalwt[complete.cases(prot.bmalwt), ]
bmalwt.mat <- subset(prot.bmalwt, select = c(-gene, -Protein.IDs))
bmalwt.mat <- t(scale(bmalwt.mat, center = TRUE, scale = TRUE))
bmalwt.pca <- prcomp(bmalwt.mat)

plot(bmalwt.pca$x[, 1], bmalwt.pca$x[, 2], pch="*")
text(bmalwt.pca$x[, 1], bmalwt.pca$x[, 2], rownames(bmalwt.pca$x))


# PCA for WT CRY ----------------------------------------------------------

jgenos <- c("WT")
jtimes <- unique(subset(prot.long, geno == jgenos)$time)
prot.sub <- subset(prot.long, geno %in% jgenos & time %in% jtimes & gene != "")
prot.sub$time <- paste("ZT", prot.sub$time, sep = "")
prot.wt <- dcast(prot.sub, formula = gene + Protein.IDs ~ time, value.var = "rel.abund")
prot.wt <- prot.wt[complete.cases(prot.wt), ]
wt.mat <- subset(prot.wt, select = c(-gene, -Protein.IDs))
wt.mat <- t(scale(wt.mat, center = TRUE, scale = TRUE))
wt.pca <- prcomp(wt.mat)

plot(wt.pca$x[, 1], wt.pca$x[, 2], pch="*")
text(wt.pca$x[, 1], wt.pca$x[, 2], rownames(wt.pca$x))


# PCA for all samples -----------------------------------------------------

prot.sub <- subset(prot.long, gene != "")
prot.sub$time <- paste("ZT", prot.sub$time, sep = "")
prot.mat <- dcast(prot.sub, formula = gene + Protein.IDs ~ time + geno, value.var = "rel.abund")
prot.mat <- prot.mat[complete.cases(prot.mat), ]
prot.mat <- subset(prot.mat, select = c(-gene, -Protein.IDs))
prot.mat <- t(scale(prot.mat, center=TRUE, scale=TRUE))
prot.mat.pca <- prcomp(prot.mat)

plot(prot.mat.pca$x[, 1], prot.mat.pca$x[, 2])
text(prot.mat.pca$x[, 1], prot.mat.pca$x[, 2], rownames(prot.mat.pca$x))


# Distribution of log2 signal in Bmal WT signals --------------------------

jtimes <- unique(subset(prot.long, geno == "BmalWT")$time)
ggplot(subset(prot.long, time %in% jtimes), aes(x = rel.abund)) + 
  geom_histogram(bins = 80) + 
  facet_grid(geno ~ time) + 
  theme_bw() + 
  geom_vline(xintercept = 0) + 
  theme(aspect.ratio = 1)

