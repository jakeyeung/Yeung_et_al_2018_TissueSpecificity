# 2015-09-28
# Use BIC but on limited set of models
# Jake Yeung

w <- 2 * pi / 24

library(dplyr)
library(parallel)
library(hash)

setwd("~/projects/tissue-specificity")


# Source ------------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/MemoryManagement.R")

# Functions ---------------------------------------------------------------

GetRhythmicSharedParams2 <- function(dat.gene, rhyth.tiss, tissue.i = 2){
  # Given dat.gene and tissues that are rhythmic, return design matrices
  # for corresponding model
  w <- 2 * pi / 24
  dat.gene$is.rhyth <- apply(dat.gene, 1, function(dat.row) as.numeric(dat.row[tissue.i] %in% rhyth.tiss))
  my_mat.shared <- model.matrix(exprs ~ 0 + tissue + experiment + I(is.rhyth*(sin(w * time))), data = dat.gene)
}

ConcatenateRhythmicTissues <- function(des.mat, des.mat.rhyth.sin, des.mat.rhyth.cos, rhyth.tiss){
  # Concatenate rhythmic tissues so they have shared parameters
  grepl.str <- paste(rhyth.tiss, collapse = "|")
  des.mat.cols <- lapply(list(des.mat.rhyth.sin, des.mat.rhyth.cos), function(d.mat){
    d.mat.sub <- d.mat[, colnames(d.mat)[grepl(grepl.str, colnames(d.mat))]]
    d.mat.sub <- as.matrix(apply(d.mat.sub, 1, sum))
  })
  des.mat.cols <- do.call(cbind, des.mat.cols)
  des.mat <- cbind(des.mat, des.mat.cols)
  return(des.mat)
}

# Load and plot on 12 conditions -------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)

# remove WFAT
dat.long <- subset(dat.long, tissue != "WFAT")


# Create test set ---------------------------------------------------------

tissues <- as.character(unique(dat.long$tissue))
# tiss.test <- c("Liver", "Kidney", "Adr", "BFAT", "Mus", "Aorta", "Heart")
# tiss.test <- tissues
tiss.test <- c("Liver", "Lung")
jgene <- "1110025L11Rik"
# jgene <- "Nr1d1"
dat.gene <- subset(dat.long, gene == jgene & tissue %in% tiss.test)
dat.gene$tissue <- factor(as.character(dat.gene$tissue), levels = tiss.test)

start <- Sys.time()
my_fits <- MakeDesMatRunFit(dat.gene, n.rhyth.max = 3, normalize.weights = FALSE, cutoff = 0.01)
print(Sys.time() - start)


# Run for ALL genes -------------------------------------------------------

# for expressed genes only
dat.exprs <- dat.long %>%
  filter(experiment == "rnaseq") %>%
  group_by(tissue, gene) %>%
  summarise(exprs.mean = mean(exprs)) %>%
  filter(exprs.mean > 5)

set.seed(0)
genes.exprs <- unique(dat.exprs$gene)
genes.exprs.sub <- sample(genes.exprs, size = 0.1 * length(genes.exprs))
# genes.exprs.sub <- genes.exprs

dat.sub <- subset(dat.long, tissue %in% tiss.test & gene %in% genes.exprs.sub)
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = tiss.test)
dat.sub$gene <- factor(as.character(dat.sub$gene), levels = genes.exprs.sub)

dat.env <- DatLongToEnvironment(dat.sub)

my_fits <- MakeDesMatRunFitEnv(dat.env, "Xrra1")

fits.all <- lapply(ls(dat.env), function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, )
}

dat.split <- split(dat.sub, dat.sub$gene)
fits.all <- lapply(dat.split, MakeDesMatRunFit, normalize.weights=TRUE, cutoff=1e-5)
# fits.all <- mclapply(dat.split, MakeDesMatRunFit, cutoff = 1e-5, mc.cores = 50)

# my_fits <- MakeDesMatRunFit(dat.split[[1]], normalize.weights = TRUE, cutoff = 1e-5)

# fits.all <- dat.sub %>%
#   group_by(gene) %>%
#   do(MakeDesMatRunFit(., normalize.weights = TRUE, cutoff = 1e-5))