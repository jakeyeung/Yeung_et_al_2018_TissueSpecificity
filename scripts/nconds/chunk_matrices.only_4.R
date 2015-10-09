# 2015-09-28
# Use BIC but on limited set of models
# Jake Yeung

set.seed(0)
w <- 2 * pi / 24

library(dplyr)
library(parallel)
library(hash)
# library(proftools)

setwd("~/projects/tissue-specificity")


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

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)

# remove WFAT
dat.long <- subset(dat.long, tissue != "WFAT")


# Create test set ---------------------------------------------------------

tissues <- as.character(unique(dat.long$tissue))
# tiss.test <- c("Liver", "Kidney", "Adr", "BFAT", "Mus", "Aorta", "Heart")
# tiss.test <- c("Liver", "Lung", "Kidney", "BFAT", "Mus")
tiss.test <- tissues


# Run for ALL genes -------------------------------------------------------

# for expressed genes only
dat.exprs <- dat.long %>%
  filter(experiment == "rnaseq") %>%
  group_by(tissue, gene) %>%
  summarise(exprs.mean = mean(exprs)) %>%
  filter(exprs.mean > 5)

genes.exprs <- unique(dat.exprs$gene)
# genes.exprs.sub <- sample(genes.exprs, size = 0.1 * length(genes.exprs))
genes.exprs.sub <- genes.exprs


# Make matrices -----------------------------------------------------------

dat.sub <- subset(dat.long, tissue %in% tiss.test & gene %in% genes.exprs.sub)
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = tiss.test)
dat.sub$gene <- factor(as.character(dat.sub$gene), levels = genes.exprs.sub)

dat.gene <- subset(dat.sub, gene == "Nr1d1")  # for making matrices
chunks <- 200000
max.n <- 4
n.tiss <- length(tiss.test)

outdir <- paste0("/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_", n.tiss, "_only_", max.n, "_chunks_", chunks)

dir.create(outdir)
start <- Sys.time()
MakeDesMatChunks(dat.gene, out.dir = outdir, tiss.test, n.rhyth.max = max.n, w = w, sparse = TRUE, chunks=chunks, only.n.params = 4)
print(Sys.time() - start)


