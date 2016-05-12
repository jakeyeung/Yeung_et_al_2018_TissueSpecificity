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

CountModels <- function(fits.best.sub){
  fits.sum <- fits.best.sub %>%
    group_by(model) %>%
    summarise(count = length(model))
  fits.sum <- fits.sum[order(fits.sum$count, decreasing = TRUE), ]
  fits.best.sub$model <- factor(as.character(fits.best.sub$model), levels = fits.sum$model)
  fits.sum <- OrderDecreasing(fits.sum, jfactor = "model", jval = "count")
  return(fits.sum)
}

PlotPolarHistogram <- function(fits.sub, countstring="Count"){
  # countstring="Count" or "NormCount"
  fits.sub.sum <- CountModels(fits.sub)
  # bin manually
  # 24 and 0 are the same bin
  fits.sub$bin <- sapply(fits.sub$phase.avg, function(p){
    if (p > 0.5){
      phase.bin <- round(p, digits = 0)
    } else {
      phase.bin <- 24
    }
    return(phase.bin)
  })
  fits.sub.forhist <- fits.sub %>%
    group_by(model, bin) %>%
    summarise(Count = length(bin))
  fits.sub.forhist$bin <- fits.sub.forhist$bin - 0.5
  
  # add label to indicate number of genes
  n.genes.hash <- hash(as.character(fits.sub.sum$model), fits.sub.sum$count)
  # rename tissues
  fits.sub.forhist$model.n <- as.factor(sapply(as.character(fits.sub.forhist$model), function(m) paste0(m, " N=", n.genes.hash[[m]])))
  tiss.ordered <- sapply(as.character(fits.sub.sum$model), function(m) paste0(m, " N=", n.genes.hash[[m]]))
  fits.sub.forhist$model.n <- factor(as.character(fits.sub.forhist$model.n), levels = tiss.ordered)
  # keep original labels
  fits.sub.forhist$tissue <- sapply(as.character(fits.sub.forhist$model.n), function(jlab){
    strsplit(jlab, " ")[[1]][[1]]
  })
  fits.sub.forhist$tissue <- factor(as.character(fits.sub.forhist$tissue), levels = names(tiss.ordered))
  fits.sub.forhist <- fits.sub.forhist %>%
    group_by(model) %>%
    mutate(NormCount = Count / sum(Count))
  
  # normalized by sum
  m2 <- ggplot(fits.sub.forhist, aes_string(x = "bin", y = countstring)) +
    geom_bar(stat = "identity", width = 1) + 
    facet_wrap(~tissue) + 
    scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) +
    expand_limits(y = 0) +
    expand_limits(x = c(0, 24)) +
    xlab("Phase (h)") + 
    ylab("Count") +
    theme_bw() + 
    theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
          axis.title = element_text(size = 24),
          strip.text = element_text(size = 12)) + 
    coord_polar(theta = "x")
  return(m2)
}

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

