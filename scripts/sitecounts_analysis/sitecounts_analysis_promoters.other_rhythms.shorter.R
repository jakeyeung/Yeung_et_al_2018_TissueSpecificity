# 2016-02-13
# Jake Yeung
# Quantify significance of sitecounts using linear discriminant analysis

rm(list=ls())


library(ggplot2)
library(reshape2)
library(dplyr)
library(hash)
library("penalizedLDA")

dist.ref <- 500  # 500 left and right of promoter is reference
# sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
# dist <- 500
sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_50000_dist_sum_multigene/sitecounts.50000.multigene.mat"
dist <- 50000  # needs rescaling
# sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_1000_dist_sum_multigene/sitecounts.1000.multigene.mat"
# dist <- 1000

# Function ----------------------------------------------------------------

source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ShannonEntropy.R")
source("scripts/functions/LdaFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/N.long.promoters_500.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/N.long.sum.bytiss.all_genes.Robj", verbose=T)

dat.mean.rnaseq <- subset(dat.long, tissue != "WFAT" & experiment == "rnaseq" & gene %in% fits.best$gene) %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))

minamp.kidliv <- 0.25
minamp.liv <- 0.5

fg.models <- "Adr"
flat.models <- ""
rhyth.models <- as.character(subset(fits.best, n.rhyth >= 8)$model)

adr.genes <- subset(fits.best, model %in% fg.models)$gene
rhyth.genes <- subset(fits.best, model %in% rhyth.models)$gene

m <- PlotMeanExprsOfModel(dat.mean.rnaseq, genes = adr.genes, jmodel = "Adr", sorted = TRUE)
m2 <- PlotOverlayTimeSeries(dat.long, genes = adr.genes, tissues = c("Adr"), jscale = TRUE, jalpha = 0.01, jtitle = "Overlay plot of Adr-rhythmic genes")
multiplot(m, m2, cols = 2)

# plot rhyth.models (should do SVD)
PlotAmpPhase(subset(fits.best, model %in% rhyth.models), amp_str = "amp.avg", phase_str = "phase.avg", lab_str = "gene")


# Rhythm vs Flat ----------------------------------------------------------

flat.mat <- LongToMat(subset(N.long, model %in% flat.models))
rhyth.mat <- LongToMat(subset(N.long, model %in% rhyth.models))

flat.v.rhyth <- DoLdaPromoters(rhyth.mat, flat.mat)
out.df <- data.frame(motif = colnames(flat.v.rhyth$x), discrim.flat = 1 * -flat.v.rhyth$discrim, discrim.rhyth = 1 * -flat.v.rhyth$discrim)

ggplot(out.df, aes(x = discrim.flat, y = discrim.rhyth, label = motif)) + geom_text() + ggtitle(paste("Rhyth model: rhythmic genes")) + 
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
  theme_bw(24) + xlab("Loadings against flat genes") + ylab("Loadings against rhythmic genes")
