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
# tiss.test <- c("Liver", "Kidney", "Adr", "BFAT", "Mus", "Aorta", "WFAT")
# tiss.test <- tissues
tiss.test <- c("Liver", "Lung")
dat.gene <- subset(dat.long, gene == "Nr1d1" & tissue %in% tiss.test)
dat.gene$tissue <- factor(as.character(dat.gene$tissue), levels = tiss.test)

my_mat <- MakeRhythmicDesignMatrices(dat.gene, simplify = TRUE)

my_fits <- FitModels(dat.gene, my_mat, get.criterion = "BIC", normalize.weights=TRUE)

mat <- my_mat[[5]]

fit.fast <- lm.fit(mat, dat.gene$exprs)
fit.lm <- lm(exprs ~ 0 + mat, dat.gene)
fit.ced <- compute_RSS(x = dat.gene$exprs, matX = mat)

bic.fast <- BICFromLmFit(fit.fast$coefficients, fit.fast$residuals)
bic.ced <- compute_BIC(fit.ced, n = length(dat.gene$exprs))
bicdiff <- BIC(fit.lm) -  BICFromLmFit(fit.fast$coefficients, fit.fast$residuals)
bicdiff

n <- length(fit.fast$residuals)
RSS <- sum(fit.fast$residuals ^ 2)
k <- length(fit.fast$coefficients)

L.lm <- print(logLik(fit.lm))
L.fast <-log(sum(fit.fast$residuals ^ 2) / n)

(criterion <- -2 * L.lm + k * log(n))  # closer
BIC(fit.lm)

# check BIC weights if differing
fits.fast <- lapply(my_mat, function(mat){
  return(lm.fit(mat, dat.gene$exprs))
})
fits.lm <- lapply(my_mat, function(mat) lm(exprs ~ 0 + mat, dat.gene))

bics.fast <- unlist(lapply(fits.fast, function(fit.f) BICFromLmFit(fit.f$coefficients, fit.f$residuals)))
bics.lm <- unlist(lapply(fits.lm, function(fit.l) BIC(fit.l)))

bics.lm.norm <- exp(-0.5 * bics.lm) / sum(exp(-0.5 * bics.lm))

bics.fast.norm <- exp(-0.5 * bics.fast) / sum(exp(-0.5 * bics.fast))

# normalize weights
weight.sum <- lapply(fit, function(f){
  f$weight
})

anova(fit[[5]], fit[[4]])

