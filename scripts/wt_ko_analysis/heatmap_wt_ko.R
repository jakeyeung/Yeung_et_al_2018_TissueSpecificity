# 2016-06-10
# Jake Yeung
# heatmap_wt_ko.R
# Summarize by heatmap

rm(list=ls())

library(dplyr)

source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

load("Robjs/wtko.dat.wtko.hog.commongenes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/wtko.fits.all.long.commongenes.Robj", v=T)

# heatmap -----------------------------------------------------------------

problem.genes <- c("Ccl25")

genes.liv <- as.character(subset(fits.best, model == "Liver")$gene)
fits.sub.liv <- subset(fits.all.long.wtkohog, gene %in% genes.liv & model != "")
genes.liv.wtliv <- as.character(subset(fits.sub.liv, model %in% c("WT,Liver", "WT;Liver", "Liver", "WT"))$gene)
jgenes <- subset(fits.best, model == "Liver" & ! gene %in% problem.genes)
jgenes <- as.character(jgenes[order(jgenes$phase.avg), ]$gene)

jgenes <- jgenes[jgenes %in% genes.liv.wtliv]

PlotGeneAcrossTissues(subset(dat.wtko.hog, gene == "Lars"))

jsub <- subset(dat.wtko.hog, gene %in% jgenes & tissue != "Liver") %>%
  group_by(gene, tissue) %>%
  mutate(exprs.scaled = scale(exprs, center = TRUE, scale = FALSE))

jsub.gene <- subset(jsub, gene == "Slc16a4")
jsub.gene$exprs <- jsub.gene$exprs.scaled
PlotGeneAcrossTissues(jsub.gene)

M <- dcast(subset(jsub, experiment != "rnaseq" & tissue != "Liver"), formula = gene ~ tissue + time, value.var = "exprs.scaled")
rownames(M) <- M$gene; M$gene <- NULL
M <- M[jgenes, ]
M <- M[complete.cases(M), ]

PlotExprsHeatmap(M, jtitle = "", blackend = 0.1, minval = -1, maxval = 1, jdendro = "none")
