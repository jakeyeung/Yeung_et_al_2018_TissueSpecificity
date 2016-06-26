# 2016-06-23
# Jake Yeung

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


# Load --------------------------------------------------------------------

# indir <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds/removesampsFALSE"
# indir <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds"
indir <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/nconds/staggered_timepoints"
outrobj <- "Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.Robj"
inf <- list.files(indir, pattern = "*.Robj")

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

# save(fits.long, file="Robjs/liver_kidney_atger_nestle/fits.long.multimethod.Robj")

# See clock genes ---------------------------------------------------------

# clockgenes <- GetClockGenes()
# fits.cg <- subset(fits.long, gene %in% clockgenes)


# Filter lowly expressed genes --------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)

genes.all <- unique(as.character(dat.long$gene))
dat.long <- RemoveLowExprsPseudoShortGenes(dat.long, ggbiotype = "protein_coding", gglength = 250, jcutoff = 1, show.plot=FALSE)

genes.keep <- unique(as.character(dat.long$gene))

print(head(genes.all))
print(head(genes.keep))
genes.cut <- genes.all[which(!genes.all %in% genes.keep)]
print(head(genes.cut))



# dat.mean <- dat.long %>%
#   group_by(gene) %>%
#   summarise(exprs.max = max(exprs))
# 
# plot(density(dat.mean$exprs.max))
# jcutoff <- 1
# abline(v=jcutoff)
# 
# genes.cut <- as.character(subset(dat.mean, exprs.max <= jcutoff)$gene)


# Get best ----------------------------------------------------------------

fits.long.filt <- subset(fits.long, !gene %in% genes.cut) %>%
  group_by(gene, method) %>%
  filter(weight.raw == min(weight.raw))

print(head(fits.long.filt))

save(fits.long.filt, file=outrobj)

