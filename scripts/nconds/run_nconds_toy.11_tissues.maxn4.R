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
# tiss.test <- tissues
tiss.test <- c("Liver", "Lung", "Kidney", "BFAT", "Mus")
# tiss.test <- tissues
# jgene <- "1110025L11Rik"
# # jgene <- "Nr1d1"
# dat.gene <- subset(dat.long, gene == jgene & tissue %in% tiss.test)
# dat.gene$tissue <- factor(as.character(dat.gene$tissue), levels = tiss.test)
# 
# start <- Sys.time()
# my_fits <- MakeDesMatRunFit(dat.gene, n.rhyth.max = 3, normalize.weights = FALSE, cutoff = 0.01)
# print(Sys.time() - start)


# Run for ALL genes -------------------------------------------------------

# for expressed genes only
dat.exprs <- dat.long %>%
  filter(experiment == "rnaseq") %>%
  group_by(tissue, gene) %>%
  summarise(exprs.mean = mean(exprs)) %>%
  filter(exprs.mean > 5)

genes.exprs <- unique(dat.exprs$gene)
# genes.exprs.sub <- sample(genes.exprs, size = 0.1 * length(genes.exprs))
# genes.exprs.sub <- genes.exprs
genes.exprs.sub <- GetClockGenes()[1:5]

dat.sub <- subset(dat.long, tissue %in% tiss.test & gene %in% genes.exprs.sub)
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = tiss.test)
dat.sub$gene <- factor(as.character(dat.sub$gene), levels = genes.exprs.sub)

print("Chunking data to environment")
start <- Sys.time()
dat.env <- DatLongToEnvironment(dat.sub)
print(Sys.time() - start)

# Rprof()
fits.all <- lapply(ls(dat.env), function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, tiss.test, n.rhyth.max = 4, w = 2 * pi / 24, criterion = "BIC", normalize.weights = TRUE, cutoff = 1e-5, top.n = 10, sparse = TRUE)
})
# Rprof(NULL); RprofData <- readProfileData("Rprof.out")
# flatProfile(RprofData, byTotal=TRUE)
# flatProfile(RprofData, byTotal=FALSE)

print("Running many cores and mclapply")
start <- Sys.time()
fits.all <- mclapply(ls(dat.env), function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, tiss.test, n.rhyth.max = 4, w = 2 * pi / 24, criterion = "BIC", normalize.weights = TRUE, cutoff = 1e-5, top.n = 10, sparse = TRUE)
}, mc.cores = 6)
print(Sys.time() - start)
fits.all[[1]]
save(fits.all, file = "/home/yeung/projects/tissue-specificity/Robjs/nconds2/fits.Robj")
