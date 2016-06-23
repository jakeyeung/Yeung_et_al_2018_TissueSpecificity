# 2016-06-23
# Jake Yeung

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetClockGenes.R")


# Load --------------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds/removesampsFALSE"
inf <- list.files(indir)

fits.long <- expandingList()
for (f in inf){
  method <- strsplit(f, split = "\\.")[[1]][[4]]  # manually change this!
  fpath <- file.path(indir, f)
  load(fpath)
  fits.all.long$method <- method
  fits.long$add(fits.all.long)
}
fits.long <- fits.long$as.list()
fits.long <- do.call(rbind, fits.long)

# See clock genes ---------------------------------------------------------

clockgenes <- GetClockGenes()


fits.cg <- subset(fits.long, gene %in% clockgenes)
