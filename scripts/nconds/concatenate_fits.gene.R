# 2015-10-11
# Jake Yeung
# Analyze fits from vitalit by gene to allow parallelization

library(dplyr)

# setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/AppendListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")

# Functions ---------------------------------------------------------------


# Load --------------------------------------------------------------------

top.n <- 100
args <- commandArgs(trailingOnly = TRUE)
fitdir <- args[[1]]
genename <- args[[2]]
outdir <- args[[3]]
outname <- paste0("fits_11_tiss_max_3.top_", top.n, ".", genename, ".Robj")
outpath <- file.path(outdir, outname)

start <- Sys.time()
fits <- LoadFitsFromModels(fitdir)
fits.top <- GetTopNModelsFromFits(fits, top.n)
fits.long.gene <- ListToLong(fits.top, genename = gene, top.n = top.n, period = 24)

save(fits.long.gene, file = outpath)
print(Sys.time() - start)
