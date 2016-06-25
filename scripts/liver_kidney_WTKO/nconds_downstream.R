# 2016-06-23
# Jake Yeung

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetClockGenes.R")


# Load --------------------------------------------------------------------

# indir <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds/removesampsFALSE"
indir <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds"
inf <- list.files(indir, pattern = "*.Robj")

fits.long <- expandingList()
for (f in inf){
  print(f)
  method <- strsplit(f, split = "\\.")[[1]][[6]]  # manually change this!
  fpath <- file.path(indir, f)
  load(fpath)
  fits.all.long$method <- method
  fits.long$add(fits.all.long)
}
fits.long <- fits.long$as.list()
fits.long <- rbind_all(fits.long)

save(fits.long, file="Robjs/liver_kidney_atger_nestle/fits.long.multimethod.Robj")

# See clock genes ---------------------------------------------------------

# clockgenes <- GetClockGenes()
# fits.cg <- subset(fits.long, gene %in% clockgenes)


# Filter lowly expressed genes --------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)

dat.mean <- dat.long %>%
  group_by(gene) %>%
  summarise(exprs.max = max(exprs))

plot(density(dat.mean$exprs.max))
jcutoff <- 1
abline(v=jcutoff)

genes.cut <- as.character(subset(dat.mean, exprs.max <= jcutoff)$gene)


# Get best ----------------------------------------------------------------

fits.long.filt <- subset(fits.long, !gene %in% genes.cut) %>%
  group_by(gene, method) %>%
  filter(weight.raw == min(weight.raw))

save(fits.long.filt, file="Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.Robj")

