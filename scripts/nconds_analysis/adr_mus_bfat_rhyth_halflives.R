# 2016-05-27
# Rhythmic genes of Adr Mus BFAT have sharp peaks WHY?

rm(list=ls())

library(wordcloud)
library(ggplot2)
library(dplyr)
library(reshape2)
source("scripts/functions/PlotGeneAcrossTissues.R")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

dat.rnaseq <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue) %>%
  summarise(exprs = mean(exprs))

# Look at Adr -------------------------------------------------------------

adr.genes <- as.character(subset(fits.best, model == "Adr")$gene)

# Get highly expressed genes in Liver -------------------------------------

dat.liv <- subset(dat.rnaseq, tissue == "Liver")
dat.adr <- subset(dat.rnaseq, tissue == "Adr" & exprs > 5)

# Plot top 100 ------------------------------------------------------------

pdf("plots/gene_expressions/top_liver_exprs.pdf")
genes <- as.character(dat.liv[order(dat.liv$exprs, decreasing = TRUE), ]$gene)

dat.sub <- subset(dat.long, gene %in% genes[1:100])
for (i in seq(1:100)){
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == genes[i])))
}
dev.off()


# Plot lowly expressed Adrenal genes --------------------------------------

pdf("plots/gene_expressions/bottom_adr_exprs.pdf")
genes <- as.character(dat.adr[order(dat.adr$exprs, decreasing = FALSE), ]$gene)

dat.sub <- subset(dat.long, gene %in% genes[1:100])
for (i in seq(1:100)){
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == genes[i])))
}
dev.off()


# PCA of Adrenal and Liver ------------------------------------------------

mat <- dcast(subset(dat.long, experiment == "rnaseq" & tissue %in% c("Mus", "BFAT")), formula = gene ~ tissue + time, value.var = "exprs")
rownames(mat) <- mat$gene; mat$gene <- NULL
mat.pca <- prcomp(mat, center = TRUE, scale. = TRUE)

screeplot(mat.pca)

pca.i <- 4
wordcloud::textplot(x = mat.pca$rotation[, pca.i], y =  mat.pca$rotation[, pca.i + 1], words = rownames(mat.pca$rotation))


# BFAT genes analysis -----------------------------------------------------

PlotGeneAcrossTissues(subset(dat.long, gene == "Mef2c"))


