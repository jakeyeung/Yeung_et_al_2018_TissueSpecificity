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

FactorToHex <- function(i, colhash){
  return(colhash[[as.character(i)]])
}

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/fits.relamp.Robj", v=T)
fits.relamp <- subset(fits.relamp, tissue != "WFAT")

# Plot heatmaps -----------------------------------------------------------

# do for liver genes
jtiss <- "Liver"
all.tiss <- as.character(unique(fits.relamp$tissue))

ts.genes <- as.character(subset(fits.best, model == jtiss)$gene)

fits.sub <- subset(fits.relamp, gene %in% ts.genes) 
# sort by phases of jtiss
fits.sub$tissue <- factor(fits.sub$tissue, levels = c(jtiss, all.tiss[which(all.tiss != jtiss)]))
fits.sub <- fits.sub %>%
  group_by(tissue) %>%
  arrange(phase)
# order by sum of amplitude
fits.sumamp <- fits.sub %>%
  group_by(tissue) %>%
  summarise(amp.sum = sum(amp)) %>%
  arrange(desc(amp.sum))

tissue.order <- as.character(as.character(fits.sumamp$tissue))
gene.order <- as.character(subset(fits.sub, tissue == jtiss)$gene)

# fits.sub$cols <- as.factor(mapply(PhaseAmpPvalToColor, fits.sub$phase, 10, 1e-6, MoreArgs = list(rotate.hr = -8)))
fits.sub$cols <- as.factor(mapply(PhaseAmpPvalToColor, fits.sub$phase, fits.sub$amp, 0, MoreArgs = list(rotate.hr = )))
# add pvalue for alpha
fits.sub$cols.alpha <- factor(mapply(function(hex, pval) AddAlphaToHexColor(hex, alpha = SaturationCurve(-log10(pval), Vmax = 1, k = 0.25, x0 = 0)),
                              as.character(fits.sub$cols), as.numeric(fits.sub$pval)))

use.alpha <- TRUE
if (use.alpha){
  fits.sub$cols.i <- seq(nrow(fits.sub))
  colhash <- hash(as.character(fits.sub$cols.i), as.character(fits.sub$cols.alpha))
} else {
  # fits.sub$cols.i <- as.numeric(fits.sub$cols)
  fits.sub$cols.i <- seq(nrow(fits.sub))
  colhash <- hash(as.character(fits.sub$cols.i), as.character(fits.sub$cols))
}

fits.mat <- dcast(fits.sub, formula = gene ~ tissue, value.var = "cols.i")
rownames(fits.mat) <- fits.mat$gene
fits.mat$gene <- NULL
# rreoder
fits.mat <- fits.mat[gene.order, ]
fits.mat <- fits.mat[, tissue.order]
fits.mat <- as.matrix(fits.mat)

image(t(fits.mat), col = sapply(sort(unlist(fits.mat)), FactorToHex, colhash), yaxt = "n", xaxt = "n", axes=FALSE,xlab="",ylab="",srt=0)
axis(3, at = seq(0, 1, length.out = ncol(fits.mat)), labels=colnames(fits.mat), srt=0, tick=FALSE)
# axis(2, at = seq(0, 1, length.out = nrow(fits.mat)), labels=rownames(fits.mat)[1:ngenes], srt=0, tick=FALSE)

# test colors
fits.sub.test <- fits.mat[, "Liver"][800:839]
plot(y = rep(1, length(fits.sub.test)), x = seq(1, length.out = length(fits.sub.test), by = 1), col=sapply(fits.sub.test, FactorToHex, colhash), cex=2, pch = 15)

fits.sub2 <- subset(fits.sub, tissue == "Liver")[800:839, ]
plot(y = rep(1, length(fits.sub2$cols.alpha)), x = seq(1, length.out = length(fits.sub2$cols), by = 1), col=as.character(fits.sub2$cols.alpha), cex=2, pch = 15)
qplot(y = rep(1, length(fits.sub2$cols)), x = seq(1, length.out = length(fits.sub2$cols), by = 1)) + geom_point(alpha = 0.25, colour=as.character(fits.sub2$cols), size = 20)
