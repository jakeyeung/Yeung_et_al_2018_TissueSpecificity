# 2016-06-08
# Jake Yeung

rm(list=ls())

library(ggplot2)
library(dplyr)
library(reshape2)
library(hash)

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/ColorFunctions.R")



# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/fits.relamp.Robj", v=T)
fits.relamp <- subset(fits.relamp, tissue != "WFAT")

# Plot heatmaps -----------------------------------------------------------

# do for tissue-specific genes

tissues <- c("Liver", "BFAT", "Mus", "Adr")

for (jtiss in tissues){
  PlotHeatmapAmpPhasePval(gene.list = as.character(subset(fits.best, model == jtiss)$gene), fits.relamp, use.alpha = FALSE, title = "", single.tissue = TRUE)
}


# Heatmap gene expression -------------------------------------------------


jtiss <- "Liver"
jtiss <- "Adr"
jtiss <- "BFAT"
jtiss <- "Mus"
ts.genes <- as.character(subset(fits.best, model == jtiss)$gene)
gene.order <- subset(fits.relamp, tissue == jtiss & gene %in% ts.genes)
gene.order <- as.character(gene.order[order(gene.order$phase), ]$gene)
PlotHeatmapGeneList(gene.order, dat.long, blackend = 1, minval = -2, maxval = 2, jtitle=jtiss)


# axis(2, at = seq(0, 1, length.out = nrow(fits.mat)), labels=rownames(fits.mat)[1:ngenes], srt=0, tick=FALSE)
# 
# # test colors
# fits.sub.test <- fits.mat[, "Liver"][800:839]
# plot(y = rep(1, length(fits.sub.test)), x = seq(1, length.out = length(fits.sub.test), by = 1), col=sapply(fits.sub.test, FactorToHex, colhash), cex=2, pch = 15)
# 
# fits.sub2 <- subset(fits.sub, tissue == "Liver")[800:839, ]
# plot(y = rep(1, length(fits.sub2$cols.alpha)), x = seq(1, length.out = length(fits.sub2$cols), by = 1), col=as.character(fits.sub2$cols.alpha), cex=2, pch = 15)
# qplot(y = rep(1, length(fits.sub2$cols)), x = seq(1, length.out = length(fits.sub2$cols), by = 1)) + geom_point(alpha = 0.25, colour=as.character(fits.sub2$cols), size = 20)
