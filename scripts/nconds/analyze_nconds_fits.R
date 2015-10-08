# 2015-10-07
# Jake Yeung
# Analyze fits from vitalit

library(dplyr)
library(ggplot2)
library(reshape2)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/AppendListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")

load("Robjs/dat.long.fixed_rik_genes.Robj")


# see individual gene -----------------------------------------------------

genedir <- "/home/yeung/projects/tissue-specificity/results/nconds/fits_11_tiss_3_max_cluster"


# Load --------------------------------------------------------------------

outpath <- "/home/yeung/projects/tissue-specificity/results/nconds/fits_11_tiss_max_3.top_50.Robj"
# outpath <- "/home/yeung/projects/tissue-specificity/results/nconds/fits_11_tiss_max_3.Robj"
outpath <- "/home/yeung/projects/tissue-specificity/results/nconds/fits_11_tiss_max_3.top_99.Robj"
load(outpath, verbose = T)

fits.long$n.params <- sapply(fits.long$model, function(m){
  length(strsplit(as.character(m), ";")[[1]])
})

ggplot(fits.long, aes(x = weight)) + geom_density() + facet_wrap(~n.rhyth)


# Get top fits from models ------------------------------------------------

fits.long.top <- fits.long %>%
  group_by(gene) %>%
  filter(weight == max(weight))

# sink("/home/yeung/projects/tissue-specificity/results/nconds/genes_with_max_3/genes_with_max_3.txt")
# for (gene in subset(fits.long.top, n.params == 3)$gene){
#   cat(gene)
#   cat("\n")
# }
# sink()

ggplot(fits.long.top, aes(x = weight)) + geom_density() + facet_wrap(~n.rhyth)


# Check my fits are proper ------------------------------------------------

tissues <- as.character(unique(dat.long$tissue))
tissues <- tissues[which(tissues != "WFAT")]
dat.gene <- subset(dat.long, gene == jgene & tissue %in% tissues)
dat.gene$tissue <- factor(as.character(dat.gene$tissue), levels = tissues)

jgene <- "Rgs16"
# jmodel <- "Liver;Adr,BFAT,Heart,Hypo"
jmodel <- "Liver"

X <- MakeDesMatFromModelName(dat.gene, model.name = jmodel, tissues)

fit.gene <- lm.fit(X, dat.gene$exprs)
fit.bic <- BICFromLmFit(fit.gene$coefficients, fit.gene$residuals)
print(fit.bic)
print(exp(-0.5 * fit.bic))
