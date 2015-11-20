# Jake Yeung
# tissue_specific_paper.supplemental.R
# What are mean expression of tissue modules?  
# 2015-11-19

library(dplyr)
library(ggplot2)
library(reshape2)
library(hash)

setwd("/home/yeung/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots_paper"

# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")


# Functions ---------------------------------------------------------------



# Load --------------------------------------------------------------------

# load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")
load("Robjs/dat.long.fixed_rik_genes.Robj")


# PCA on mean expression --------------------------------------------------

dat.mean <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))


# Count of genes per model ------------------------------------------------

fits.count <- fits.best %>%
  group_by(model) %>%
  summarise(count = length(gene))
fits.count <- fits.count[order(fits.count$count, decreasing = T), ]


# Take all adrenal-rhythmic genes and plot their rhythms ------------------

adr.genes <- as.character(subset(fits.best, model == "Adr")$gene)

m <- PlotOverlayTimeSeries(dat.long, adr.genes, tissues = "Adr", jalpha = 0.05, jtitle = "Adrenal-specific rhythmic genes")
print(m)

# BFAT --------------------------------------------------------------------

BFAT.genes <- as.character(subset(fits.best, model == "BFAT")$gene)

m <- PlotOverlayTimeSeries(dat.long, BFAT.genes, tissues = "BFAT", jalpha = 0.05, jtitle = "BFAT-specific rhythmic genes")
print(m)


# Mus ---------------------------------------------------------------------

jtiss <- "Mus"
Mus.genes <- as.character(subset(fits.best, model == jtiss)$gene)

m <- PlotOverlayTimeSeries(dat.long, Mus.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)



# Liver -------------------------------------------------------------------


jtiss <- "Liver"
Mus.genes <- as.character(subset(fits.best, model == jtiss)$gene)

m <- PlotOverlayTimeSeries(dat.long, Mus.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)


# Is Liver-specific genes liver-specific by mean exprs? -------------------

# jtiss <- "Liver"
# Liver.genes <- as.character(subset(fits.best, model == jtiss)$gene)
# 
# dat.mean.sub <- subset(dat.mean, gene %in% Liver.genes)
# dat.mat <- dcast(dat.mean.sub, gene ~ tissue, value.var = "exprs.mean")
# rownames(dat.mat) <- dat.mat$gene; dat.mat$gene <- NULL
# 
# dat.pca <- prcomp(t(dat.mat), center = T, scale. = T)
# screeplot(dat.pca, type = "lines")
# barplot(dat.pca$x[, 1])
# 
# # boxplot of genes
# ggplot(dat.mean.sub, aes(x = tissue, y = exprs.mean)) + geom_boxplot() + theme_bw(24) + 
#   theme(aspect.ratio=1,
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   xlab("") + ylab("Mean expression of gene") +
#   ggtitle(paste0("Mean expression level of genes in ", jtiss, " module"))


# What about Adr BFAT Mus -------------------------------------------------

jtiss <- "Liver"
Liver.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, Adr.genes, jtiss)

jtiss <- "Adr"
Adr.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, Adr.genes, jtiss)

jtiss <- "BFAT"
BFAT.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, BFAT.genes, jtiss)

jtiss <- "Mus"
Mus.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, Mus.genes, jtiss)

jtiss <- "Adr;Aorta;BFAT"
BFAT.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, BFAT.genes, jtiss)

jtiss <- "Aorta;BFAT"
BFAT.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, BFAT.genes, jtiss)

jtiss <- "Aorta;BFAT;Mus"
BFAT.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean, BFAT.genes, jtiss)


