# 2015-10-30
# Jake Yeung
# Test hypothesis that BIC weights are convex
library(dplyr)
library(ggplot2)
library(reshape2)


# Functions ---------------------------------------------------------------

NMin <- function(vec){
  # check number of minimums in an ordered vector
  # algorithm: check left and right neighbhours. If both larger, then we are at minimum
  n.min <- 0
  for (i in 2:(length(vec) - 1)){
    if (vec[i] < vec[i - 1] & vec[i] < vec[i + 1]){
      n.min <- n.min + 1
    }
  }
  return(n.min)
}

# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/OuterComplex.R")

omega <- 2 * pi / 24
filt.tiss <- c("WFAT")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/dat.complex.fixed_rik_genes.Robj")

load("/home/yeung/projects/tissue-specificity/Robjs/fits.best.5_tiss.nparams.Robj", verbose=T)  # fits.best.nparam

ggplot(subset(fits.best.nparam, gene == "Arntl"), aes(x = n.params, y = weight.raw)) + geom_line()


# Summarise number of minimum ---------------------------------------------

fits.nmin <- fits.best.nparam %>%
  group_by(gene) %>%
  summarise(nmin = NMin(weight.raw))

badgenes <- subset(fits.nmin, nmin > 1)$gene
ggplot(subset(fits.best.nparam, gene %in% badgenes), aes(x = n.params, y = weight.raw)) + geom_line() + facet_wrap(~gene)


