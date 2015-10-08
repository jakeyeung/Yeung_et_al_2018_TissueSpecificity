# Load chunk, datgene, then run nconds
# 2015-10-06

set.seed(0)
w <- 2 * pi / 24

library(dplyr)
library(parallel)
library(hash)

setwd("~/projects/tissue-specificity")

writedir <- "/home/yeung/projects/tissue-specificity/data/nconds2/datlong_11_tiss_by_gene"
# Source ------------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")
source("scripts/functions/GetClockGenes.R")

# Load and plot on 12 conditions -------------------------------------------

# BEGIN: MAKE DAT GENE CHUNKS
# load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)
# 
# # remove WFAT
# dat.long <- subset(dat.long, tissue != "WFAT")
# 
# # Create test set ---------------------------------------------------------
# 
# tissues <- as.character(unique(dat.long$tissue))
# tiss.test <- tissues
# 
# 
# # Get genes ---------------------------------------------------------------
# 
# genes.exprs.sub <- GetClockGenes()
# 
# dat.sub <- subset(dat.long, tissue %in% tiss.test & gene %in% genes.exprs.sub)
# dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = tiss.test)
# dat.sub$gene <- factor(as.character(dat.sub$gene), levels = genes.exprs.sub)
# 
# ChunkDatGenesToFile(dat.sub, writedir)
# END: MAKE DAT GENE CHUNKS 

# Load datgene and matrix chunk, run my fits ------------------------------

chunkpath <- "/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_max4_overnight/chunk.2039.Robj"
chunkpath <- "/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_max4_chunks_of_100000/chunk.1.Robj"
chunkpath <- "/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_max_3_chunks_2e+05/chunk.1.Robj"
load(chunkpath, verbose=T)

datgenepath <- "/home/yeung/projects/tissue-specificity/data/nconds2/datlong_11_tiss_by_gene/Arntl.Robj"
load(datgenepath, verbose=T)

start <- Sys.time()
fits <- LoadDesMatDatGeneRunFits(dat.gene, des.mats.list, criterion = "BIC", normalize.weights = TRUE, top.n = 10, sparse = TRUE)
print(Sys.time() - start)


# Fit all -----------------------------------------------------------------

chunkdir <- "/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_max4_overnight"
outdir <- "/home/yeung/projects/tissue-specificity/data/nconds2/fits_11_tiss_chunks"

lapply(list.files(chunkdir), function(chunkpath){
  load(datgenepath, verbose=F)
  chunk.id <- strsplit(chunkpath, "\\.")[[1]][[2]]
  fits <- LoadDesMatDatGeneRunFits(dat.gene, des.mats.list, criterion = "BIC", normalize.weights = TRUE, top.n = 10, sparse = TRUE)
  save(fits, file = file.path(outdir, paste0("chunks.", chunk.id, ".fit.Robj")))
})

X <- as.matrix(des.mats.list[[202]]$mat)
Y <- dat.gene$exprs
start <- Sys.time()
for (i in seq(10000)){
  fit <- lm.fit(X, Y)
}
print(Sys.time() - start)

