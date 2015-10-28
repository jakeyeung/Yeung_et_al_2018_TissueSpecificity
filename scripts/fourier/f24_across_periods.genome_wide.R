# f24_across_periods.R
# Jake Yeung
# 2015-09-22
# library("devtools")
# dev_mode()
# 
# install("~/projects/f24")  # use jake branch
# library(f24.R2.cycling)

library(dplyr)
library(ggplot2)
setwd("/home/yeung/projects/tissue-specificity")

args <- commandArgs(trailingOnly=TRUE)
p.min <- as.numeric(args[1])
p.max <- as.numeric(args[2])
fout <- paste0("dat.fit.scan_periods.genome_wide.", p.min, "_", p.max, ".Robj")
fdir <- "Robjs"
foutpath <- file.path(fdir, fout)
print(paste("Saving to:", foutpath))

# Source ------------------------------------------------------------------

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FitRhythmicAcrossPeriods.R")

# Load data ---------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
# load("Robjs/dat.long.Robj")


# Perform f24 -------------------------------------------------------------

periods <- seq(p.min, p.max, by = 0.1)
# across all genes
start <- Sys.time()
dat.fit.periods.genome_wide <- FitRhythmicScanPeriods(dat.long, periods, cores = 50)
print(Sys.time() - start)
head(dat.fit.periods.genome_wide)
save(dat.fit.periods.genome_wide, file = foutpath)
# 2.5 hours later...
