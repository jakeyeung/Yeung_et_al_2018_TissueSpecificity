# Jake Yeung
# Date of Creation: 2018-01-20
# File: ~/projects/tissue-specificity/scripts/kidney_WT_KO/downstream_nconds.R
# Downstream analysis across many g's

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/BiomartFunctions.R")


# Constants ---------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds_kidney_WTKO"
inf <- list.files(indir, pattern = "*.Robj")
outrobj <- "Robjs/liver_kidney_atger_nestle/nconds_kidney_WTKO/summary/fits.kidney_WTKO.multimethod.long.filtbest.staggeredtimepts.bugfixed.Robj"

# Load --------------------------------------------------------------------


fits.long <- expandingList()
for (f in inf){
  print(f)
  method <- strsplit(f, split = "\\.")[[1]][[8]]  # manually change this!
  print(method)
  fpath <- file.path(indir, f)
  load(fpath)
  fits.all.long$method <- method
  fits.long$add(fits.all.long)
}
fits.long <- fits.long$as.list()
fits.long <- rbind_all(fits.long)


# Filter lowly expressed genes --------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)

dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)

dat.long <- subset(dat.long, tissue == "Kidney_SV129" | tissue == "Kidney_BmalKO")
dat.long$tissue <- factor(as.character(dat.long$tissue), levels = c("Kidney_SV129", "Kidney_BmalKO"))

genes.all <- unique(as.character(dat.long$gene))
dat.long <- RemoveLowExprsPseudoShortGenes(dat.long, ggbiotype = "protein_coding", gglength = 0, jcutoff = 1, show.plot=FALSE)

genes.keep <- unique(as.character(dat.long$gene))

print(head(genes.all))
print(head(genes.keep))
genes.cut <- genes.all[which(!genes.all %in% genes.keep)]
print(head(genes.cut))

# Get best ----------------------------------------------------------------

fits.long.filt <- subset(fits.long, !gene %in% genes.cut) %>%
  group_by(gene, method) %>%
  filter(weight.raw == min(weight.raw))

print(head(fits.long.filt))

save(fits.long.filt, file=outrobj)

