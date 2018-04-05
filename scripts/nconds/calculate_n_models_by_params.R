# Analyze rhythmic models
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
load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)

# remove WFAT
dat.long <- subset(dat.long, tissue != "WFAT")

# Create test set ---------------------------------------------------------

tissues <- as.character(unique(dat.long$tissue))
tiss.test <- c("Liver", "Kidney", "BFAT", "Mus", "Hypo", "Adr", "BS", "Aorta")

genes.exprs.sub <- GetClockGenes()
dat.sub <- subset(dat.long, tissue %in% tiss.test & gene %in% genes.exprs.sub)
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = tiss.test)
dat.sub$gene <- factor(as.character(dat.sub$gene), levels = genes.exprs.sub)

start <- Sys.time()
dat.gene <- subset(dat.sub, gene == "Nr1d1")
MakeDesMatChunks(dat.gene, out.dir = "/home/yeung/projects/tissue-specificity/data/nconds2/mats_chunks_test", tissues = tiss.test, n.rhyth.max = length(tiss.test), w = w, sparse = TRUE, 
                 chunks = 1000)
print(Sys.time() - start)


# How fast to generate matrices? ------------------------------------------
