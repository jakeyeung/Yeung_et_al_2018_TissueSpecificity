# 2016-06-27
# Jake Yeung

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(reshape2)

source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/BiomartFunctions.R")

args <- commandArgs(trailingOnly = TRUE)
jmodel <- args[[1]]
outdir <- args[[2]]
jmethod <- "g=1001"

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.Robj", v=T)

dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)

genes <- as.character(subset(fits.long.filt, model == jmodel & method == jmethod)$gene)

if (length(genes) == 0) stop("Error: Empty gene list!")

# Make mat ----------------------------------------------------------------

M <- dcast(subset(dat.long, gene %in% genes), formula = gene ~ tissue + time, value.var = "exprs")
rownames(M) <- M$gene; M$gene <- NULL

# Center rows
M.center <- sweep(M, MARGIN = 1, rowMeans(M), FUN = "-")

# Write matrices by tissue ------------------------------------------------

tissues <- as.character(unique(dat.long$tissue))

# outdir <- "/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney"
dir.create(outdir)
for (tiss in tissues){
  print(tiss)
  M.sub <- cbind(Gene.ID = rownames(M.center), M.center[, grepl(pattern = tiss, x = colnames(M.center)), ])
  write.table(M.sub, file = file.path(outdir, paste0(tiss, ".mat")), quote = FALSE, sep = "\t", row.names = FALSE)
}
