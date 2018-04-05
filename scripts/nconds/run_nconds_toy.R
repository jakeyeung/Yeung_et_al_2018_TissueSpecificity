# 2015-09-28
# Use BIC but on limited set of models
# Jake Yeung

set.seed(0)
w <- 2 * pi / 24

library(dplyr)
library(parallel)
library(hash)
library(proftools)
library(Matrix)
library(ggplot2)

setwd("~/projects/tissue-specificity")


# Source ------------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")


# Load and plot on 12 conditions -------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)

# remove WFAT
dat.long <- subset(dat.long, tissue != "WFAT")


# Create test set ---------------------------------------------------------

jtop.n <- 203
tissues <- as.character(unique(dat.long$tissue))
# tiss.test <- c("Liver", "Kidney", "Adr", "BFAT", "Mus", "Aorta", "Heart")
# tiss.test <- tissues
# tiss.test <- c("Liver", "Lung", "Kidney", "BFAT", "Mus")
tiss.test <- c("Liver", "Lung", "Kidney", "BFAT", "Heart")
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
# genes.exprs.sub <- GetClockGenes()
# genes.exprs.sub <- sample(genes.exprs, size=2000)
genes.exprs.sub <- genes.exprs

dat.sub <- subset(dat.long, tissue %in% tiss.test & gene %in% genes.exprs.sub)
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = tiss.test)
dat.sub$gene <- factor(as.character(dat.sub$gene), levels = genes.exprs.sub)


# Using environments ------------------------------------------------------


print("Chunking data to environment")
dat.env <- DatLongToEnvironment(dat.sub)

start <- Sys.time()
# Rprof()
fits.all <- mclapply(ls(dat.env), function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, tiss.test, n.rhyth.max = length(tiss.test), w = 2 * pi / 24, criterion = "BIC", normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = TRUE)
}, mc.cores = 50)
print(Sys.time() - start)

fits.all.long <- lapply(fits.all, function(x){
  gene <- x$gene
  x$gene <- NULL
  fits.temp.long <- ListToLong(x, gene, top.n = jtop.n)
})
fits.all.long <- do.call(rbind, fits.all.long)
print(Sys.time() - start)
save(fits.all.long, file = "Robjs/fits.long.5_tiss.nparams.Robj")

fits.all.long$n.params <- sapply(fits.all.long$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.all.long$n.rhyth <- sapply(fits.all.long$model, GetNrhythFromModel)
fits.all.long$amp.avg <- mapply(GetAvgAmpFromParams, fits.all.long$param.list, fits.all.long$model)
fits.all.long$phase.sd <- mapply(GetSdPhaseFromParams, fits.all.long$param.list, fits.all.long$model)
fits.all.long$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.all.long$param.list, fits.all.long$model)


# get best model from each n.params
fits.best.nparam <- fits.all.long %>%
  group_by(gene) %>%
  filter(weight.raw == min(weight.raw))

save(fits.best.nparam, file = "Robjs/fits.best.5_tiss.nparams.Robj")

# Plot BIC weight by n.params ---------------------------------------------

# ggplot(fits.best.nparam, aes(x = n.params, y = weight.raw)) + geom_line() + facet_wrap(~gene)

