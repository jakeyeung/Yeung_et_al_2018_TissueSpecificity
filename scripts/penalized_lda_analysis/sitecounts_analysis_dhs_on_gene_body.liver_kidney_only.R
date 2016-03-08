# 2016-02-19
# do penalized LDA using DHS on gene body
# redo using a new set of sitecounts matrix derived from get_tissue_spec_peaks.R

rm(list=ls())

source("scripts/functions/LdaFunctions.R")
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/dat.fit.Robj", v=T)
source("scripts/functions/PlotGeneAcrossTissues.R")

load("Robjs/N.long.sum.bytiss.all_genes.Robj", verbose=T)
# load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/fits.liver_kidney.collapsed.amp015.Robj", v=T)
fits.best <- fits.best.livkid; rm(fits.best.livkid)

library(reshape2)
library(penalizedLDA)

# Functions ---------------------------------------------------------------


# Load --------------------------------------------------------------------

minamp <- 0.25

# define parameters
# use.cross <- FALSE
sitecount.name <- "sitecount.max"
# bg.genes.type <- "Rhythmic"  # "Flat", "Rhythmic", "same" 
# bg.genes.type <- "Flat"  # "Flat", "Rhythmic", "same" 
bg.genes.type <- "same"


# jmodels <- c("Liver,Kidney")
# jmodels <- c("Liver")
jmodels <- c("Kidney")

# define foreground and background tissues
fg.tiss <- c("Kidney")
bg.tiss <- c("Liver")
# bg.tiss <- fg.tiss
print(fg.tiss)
print(bg.tiss)

out <- RunPenalizedLDA(N.long.sum.bytiss, fits.best, jmodels, fg.tiss, bg.tiss, bg.genes.type, sitecount.name, minamp = 0.25, n.bg = "all")

PlotLdaOut(out, jtitle = paste("Model:", paste(jmodels, collapse="|"), "\nFG tiss:", paste(fg.tiss, collapse="|"),"BG tiss:", paste(bg.tiss, collapse="|")))

