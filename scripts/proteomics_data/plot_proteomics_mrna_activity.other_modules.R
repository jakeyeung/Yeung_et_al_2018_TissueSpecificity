# 2016-08-16
# Jake Yeung
# explore_proteomics_data.R

rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)

source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/ProteomicsFunctions.R")

# Functions ---------------------------------------------------------------

prot.long <- LoadProteomicsData()


# Get activities ----------------------------------------------------------


load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == "g=1001")
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)
fits.long.filt$amp.avg <- mapply(GetAvgAmpFromParams, fits.long.filt$param.list, fits.long.filt$model)

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)
dat.orig <- dat.long
dat.long <- StaggeredTimepointsLivKid(dat.long)

jmod <- "Liver_BmalKO"
jtiss <- "Liver"

gene.prot <- ""
gene.lst <- list()
gene.dat <- "Nfya"
gene.act <- "NFY.A.B.C"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act)

gene.dat <- "Nfyc"
gene.act <- "NFY.A.B.C"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act)

gene.dat <- "Onecut1"
gene.act <- "ONECUT1.2"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act)

gene.dat <- "Onecut2"
gene.act <- "ONECUT1.2"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act)

# gene.dat <- "Irf1"
# gene.act <- "IRF1.2.7"
# gene.dat <- "Nfya"
# gene.act <- "NFY.A.B.C."
# gene.prot <- ""

outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
outmain <- file.path(outbase, paste0("promoters.", jmod, ".g=1001"))
indir <- file.path(outmain, "atger_with_kidney.bugfixed")
source("scripts/functions/LoadActivitiesLong.R")
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)

for (genes in gene.lst){
  gene.dat <- genes[1]; gene.act <- genes[2]
  m <- PlotmRNAActivityProtein(dat.long, act.long, prot.long = NA, gene.dat = gene.dat, gene.act = gene.act, gene.prot = gene.prot, jtiss = jtiss); print(m)
  m <- PlotmRNAActivityProtein(dat.long, act.long, prot.long = NA, gene.dat = gene.dat, gene.act = gene.act, gene.prot = gene.prot, jtiss = "both"); print(m + theme(strip.text = element_blank()))
}



