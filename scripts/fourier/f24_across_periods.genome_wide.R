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



# Source ------------------------------------------------------------------

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FitRhythmicAcrossPeriods.R")

# Load data ---------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
# load("Robjs/dat.long.Robj")


# Perform f24 -------------------------------------------------------------

# on a single gene
jgene <- "Nr1d1"
jtissue <- "Liver"
clockgenes <- GetClockGenes()

period <- 24
dat.fit <- subset(dat.long, gene %in% clockgenes) %>%
  group_by(gene, tissue) %>%
  do(FitRhythmic(., T = period, get.residuals = TRUE))
dat.fit
# dat.fit <- FitRhythmic(dat.long)

p.min <- 10
p.max <- 30
periods <- seq(p.min, p.max, by = 0.1)
# across all genes
start <- Sys.time()
dat.fit.periods.genome_wide <- FitRhythmicScanPeriods(dat.long, periods, cores = 50)
print(Sys.time() - start)
head(dat.fit.periods.genome_wide)
save(dat.fit.periods.genome_wide, file = "Robjs/dat.fit.scan_periods.genome_wide.10_to_30.Robj")
# 2.5 hours later...
