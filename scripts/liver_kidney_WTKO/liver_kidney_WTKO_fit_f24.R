# 2016-06-22
# Jake Yeung
# liver_kidney_WTKO_explore.R

rm(list=ls())
setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/FitRhythmic.R")

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)


# Fit genes ---------------------------------------------------------------

dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)

# PlotGeneAcrossTissues(subset(dat.long, gene == "Dbp"))


# Fit rhythmic parameters -------------------------------------------------

fits.bytiss <- FitRhythmicDatLong(dat.long)
save(fits.bytiss, file="Robjs/liver_kidney_atger_nestle/fits.bytiss.Robj")