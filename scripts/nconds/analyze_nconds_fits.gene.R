# 2015-10-08
# Jake Yeung
# Analyze fits gene by gene

library(dplyr)
library(ggplot2)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/AppendListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")

# load("Robjs/dat.long.fixed_rik_genes.Robj")


# Load --------------------------------------------------------------------

jgene <- "Slc34a2"
jgene <- "Dbp"
jgene <- "Rgs16"
jgene <- "Arntl"
jgene <- "Per3"
jgene <- "Nr1d1"
jgene <- "Myh7"

dirmain <- "/home/yeung/projects/tissue-specificity/data/nconds2/fits_11_tiss_chunks.max_3_test"
top.n <- 100

fitdir <- file.path(dirmain, jgene)
fits <- LoadFitsFromModels(fitdir)
fits.top <- GetTopNModelsFromFits(fits, top.n)
fits.long.gene <- ListToLong(fits.top, genename = jgene, top.n = top.n, period = 24)


# Calculate n.params and n.tissue -----------------------------------------

fits.long.gene$n.params <- sapply(fits.long.gene$model, function(m) length(strsplit(as.character(m), ";")[[1]]))

best.n <- fits.long.gene %>%
  group_by(n.params) %>%
  filter(weight == max(weight))


# Run concatenate fits and analyze ----------------------------------------

load("/home/yeung/projects/tissue-specificity/results/nconds/fits_11_tiss_clock_genes/fits_11_tiss_max_3.top_100.Robj", verbose = T)

fits.long$n.param <- sapply(fits.long$model, function(model) length(strsplit(as.character(model), ";")[[1]]))







