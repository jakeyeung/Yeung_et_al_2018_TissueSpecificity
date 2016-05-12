# 2015-11-24
# plot_tissue_specific_modules.R

rm(list=ls())
filt.tiss <- c("WFAT")

library(ggplot2)
library(dplyr)
library(hash)
library(gplots)

source("scripts/functions/PlotFunctions.R")
source("scripts/functions/LongToMat.R")
source("scripts/functions/SvdFunctions.R")


# Functions to make life easy ---------------------------------------------



# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")
load("Robjs/dat.complex.fixed_rik_genes.Robj")


# Plot histograms ---------------------------------------------------------

fits.tspec.sum <- CountModels(subset(fits.best, n.rhyth == 1)) 
m1 <- ggplot(fits.tspec.sum, aes(x = model, y = count)) + geom_bar(stat = "identity") + 
  xlab("") + 
  ylab("Count") + 
  theme_bw(24) + 
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
jmods <- c("Liver", "Adr", "BFAT", "Mus")
m2 <- PlotPolarHistogram(subset(fits.best, n.rhyth == 1 & model %in% jmods), "Count")
multiplot(m1, m2, cols = 2) 

# Plot SVD -------------------------------------------------------------------

# tissues <- unique(dat.long$tissue)
# 
# for (tiss in tissues){
#   fits.ts <- subset(fits.best, n.rhyth == 1 & model == tiss)
#   genes.ts <- as.character(fits.ts$gene)
#   print(paste("N genes for", tiss, ":", length(genes.ts)))
#   if (length(genes.ts) <= length(tissues)){
#     # problems with too few genes skip
#     next
#   }
#   
#   s.ts <- SvdOnComplex(subset(dat.complex, gene %in% genes.ts), value.var = "exprs.transformed")
#   eigens.ts <- GetEigens(s.ts, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
#   jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
#   multiplot(eigens.ts$u.plot, eigens.ts$v.plot, layout = jlayout)    
# }

