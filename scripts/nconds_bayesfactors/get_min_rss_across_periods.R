# 2016-06-12
# Jake Yeung 

library(dplyr)
library(ggplot2)
setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FitRhythmicAcrossPeriods.R")
source("scripts/functions/LiverKidneyFunctions.R")

# Load --------------------------------------------------------------------

dat <- LoadLivKid()

# Run across periods ------------------------------------------------------

p.min <- 9
p.max <- 31
dir.create("Robjs/liver_kidney")
foutpath <- "Robjs/liver_kidney/dat.fit.periods.genome_wide.Robj"

periods <- seq(p.min, p.max, by = 0.5)
# across all genes
start <- Sys.time()
dat.fit.periods.genome_wide <- FitRhythmicScanPeriods(dat, periods, cores = 50)
print(Sys.time() - start)
head(dat.fit.periods.genome_wide)
save(dat.fit.periods.genome_wide, file = foutpath)
# 2.5 hours later...