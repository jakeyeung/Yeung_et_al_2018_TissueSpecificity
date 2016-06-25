# 2016-06-25
# Jake Yeung
# First step to getting motifs from gene list 

rm(list=ls())

library(dplyr)
library(ggplot2)


outdir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/liver_kidney_wtko_modules"

# Sources -----------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")


# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.freq.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)


# Get gene lists ----------------------------------------------------------


models <- c("Kidney_SV129,Liver_SV129", "Kidney_SV129", "Liver_SV129")

for (m in models){
  outname <- paste0(gsub(",", "-", m), ".list")
  outf <- file.path(outdir, outname)
  genes <- as.character(subset(fits.long.filt, model == m)$gene)
  sink(file = outf)
  for (g in genes){
    cat(g)
    cat("\n")
  }
  sink()
}

jgene <- "Zfp651"
PlotGeneTissuesWTKO(subset(dat.long, gene == jgene))
