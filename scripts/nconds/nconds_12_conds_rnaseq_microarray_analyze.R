# 2015-09-29
# Jake Yeung

library(dplyr)
library(ggplot2)
library(reshape2)

# Functions ---------------------------------------------------------------

GetTissuesFromCoefFit <- function(fit.coef.names, tissue.colname = "tissue"){
  # Get tissues involved by grepping everything after tissue.colname
  match.names.i <- grepl(tissue.colname, fit.coef.names) & !grepl(":", fit.coef.names)
  match.names <- fit.coef.names[match.names.i]
  tissues <- sapply(match.names, function(jname) strsplit(jname, split = tissue.colname)[[1]][[2]], USE.NAMES = FALSE)
  return(tissues)
}

ExtraParamsFromFit <- function(fit.coef, tissues){
  # From coefficients of fit, extract parameters
  # colnames should contain tissue, RNA-Seq intercept, phase, amplitude, and BIC
  cnames <- tissues
  # TODO
}

SubsetByMaxBicWeight <- function(dat){
  max.i <- which.max(dat$bicweight)
  return(dat[max.i, ])
}

# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("/home/yeung/projects/tissue-specificity/Robjs/nconds.rnaseq_microarray.independent.bckup.Robj", verbose = T)

tissues <- GetTissuesFromCoefFit(names(fits.all$mod[[1]]$fit))

# check for bad fits
# fits.sum <- fits.all %>%
#   group_by(gene) %>%
#   do(SubsetByMaxBicWeight(.))
# ggplot(fits.sum, aes(x = MaxBicW)) + geom_density()

bicmat <- dcast(fits.all, formula = gene ~ rhyth.tiss, fun.aggregate = max, fill = 0, value.var = "bicweight")
rownames(bicmat) <- bicmat$gene; bicmat$gene <- NULL

# remove rows where Var.2 (flat model) is > 0
bicmat <- bicmat[which(bicmat$Var.2 == 0), ]; bicmat$Var.2 <- NULL

PlotRelampHeatmap(bicmat, jtitle = "heatmap", blackend = 0.5, yellowstart = 0.51, maxval = 1, dist.method = "manhattan")


# Get clusters ------------------------------------------------------------

clusters <- kmeans(bicmat, centers = 5)

