# 2016-06-27
# Jake Yeung

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/BiomartFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)


# Filter genes and combine tissue and genotype ----------------------------

dat.long <- RemoveLowExprsPseudoShortGenes(dat.long)
dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)

# Make mat ----------------------------------------------------------------

M <- dcast(subset(dat.long, !is.na(gene)), formula = gene ~ tissue + time, value.var = "exprs")
rownames(M) <- M$gene; M$gene <- NULL

# Center rows
M.center <- sweep(M, MARGIN = 1, rowMeans(M), FUN = "-")


# Write matrices by tissue ------------------------------------------------

tissues <- as.character(unique(dat.long$tissue))

outdir <- "/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney"
for (tiss in tissues){
  M.sub <- cbind(Gene.ID = rownames(M.center), M.center[, grepl(pattern = tiss, x = colnames(M.center)), ])
  write.table(M.sub, file = file.path(outdir, paste0(tiss, ".mat")), quote = FALSE, sep = "\t", row.names = FALSE)
}
