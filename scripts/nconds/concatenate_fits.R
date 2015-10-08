# 2015-10-07
# Jake Yeung
# Analyze fits from vitalit

library(dplyr)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/AppendListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")

load("Robjs/dat.long.fixed_rik_genes.Robj")

# Functions ---------------------------------------------------------------


# Load --------------------------------------------------------------------

top.n <- 99
# genedir <- "/home/yeung/projects/tissue-specificity/data/nconds2/fits_11_tiss_chunks.vitalit"
genedir <- "/home/yeung/projects/tissue-specificity/results/nconds/fits_11_tiss_3_max_cluster"
outdir <- "/home/yeung/projects/tissue-specificity/results/nconds"
outname <- paste0("fits_11_tiss_max_3.top_", top.n, ".Robj")
outpath <- file.path(outdir, outname)

start <- Sys.time()
fits.long.list <- expandingList()
for (gene in list.files(genedir)){
  # gene <- "9030624J02Rik"
  # gene <- "0610009L18Rik"
  gene <- "Rgs16"
  # print(gene)
  fitdir <- file.path(genedir, gene)
  fits <- LoadFitsFromModels(fitdir)
  fits.top <- GetTopNModelsFromFits(fits, top.n)
  fits.long.gene <- ListToLong(fits.top, genename = gene, top.n = top.n, period = 24)
  fits.long.list$add(fits.long.gene)
  break
}
fits.long.list <- fits.long.list$as.list()
fits.long <- do.call(rbind, fits.long.list)
head(fits.long)

fits.long$n.rhyth <- sapply(fits.long$model, GetNrhythFromModel)

# save(fits.long, file = outpath)
print(Sys.time() - start)
