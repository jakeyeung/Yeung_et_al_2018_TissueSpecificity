# 2016-02-19
# do penalized LDA using DHS on gene body
# redo using a new set of sitecounts matrix derived from get_tissue_spec_peaks.R

rm(list=ls())

source("scripts/functions/LdaFunctions.R")
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/dat.fit.Robj", v=T)
source("scripts/functions/PlotGeneAcrossTissues.R")

load("Robjs/N.long.sum.bytiss.all_genes.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
# load("Robjs/fits.liver_kidney.collapsed.amp015.Robj", v=T)

library(reshape2)
library(penalizedLDA)

# Functions ---------------------------------------------------------------


# Load --------------------------------------------------------------------

minamp <- 0.25

# define parameters
# use.cross <- FALSE
sitecount.name <- "sitecount.max"
# bg.genes.type <- "Rhythmic"  # "Flat", "Rhythmic", "same" 
bg.genes.type <- "Flat"  # "Flat", "Rhythmic", "same" 

# define rhythmic model
# jmodels <- c("Kidney")

# Kidney Liver models
# jgrep <- "^Kidney;Liver$|^Kidney,Liver$"
# jmodels <- unique(fits.best[grepl(jgrep, fits.best$model), ]$model)

jmodels <- c("Liver")
# jmodels <- c("Adr;Liver")

# define foreground and background tissues
all.tiss <- c('Cere', 'Heart', 'Kidney', 'Liver', 'Lung', 'Mus')
fg.tiss <- c("Liver")
# fg.tiss <- c("Liver")
# fg.tiss <- c("Kidney")
# fg.tiss <- c("Heart")
# bg.tiss <- all.tiss[which(! all.tiss %in% fg.tiss)]
# bg.tiss <- all.tiss[which(! all.tiss %in% fg.tiss & all.tiss != "Kidney")]
bg.tiss <- c("Kidney")
# bg.tiss <- c("Liver")
print(fg.tiss)
print(bg.tiss)

out <- RunPenalizedLDA(N.long.sum.bytiss, fits.best, jmodels, fg.tiss, fg.tiss, bg.genes.type, sitecount.name, minamp = 0.25, n.bg = "all")

PlotLdaOut(out)

