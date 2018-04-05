# 2016-08-06
# Jake Yeung
# Find ROR and HNF4 sitecounts in promoters of liver rhythmic genes

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(hash)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/GetTFs.R")



# Load --------------------------------------------------------------------


# load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
dat.orig <- dat.long

dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)
# dat.long <- SameTimepointsLivKid(dat.long)

# filter NA changes
dat.long <- subset(dat.long, !is.na(gene))


# Filter to common genes --------------------------------------------------

genes.keep <- unique(as.character(fits.long.filt$gene))

dat.long <- subset(dat.long, gene %in% genes.keep)

# Project to Frequency ----------------------------------------------------


omega <- 2 * pi / 24
dat.freq <- dat.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))

s <- SvdOnComplex(dat.freq, value.var = "exprs.transformed")

for (i in seq(1)){
  eigens <- GetEigens(s, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
}

# All periods -------------------------------------------------------------

periods <- rep(48, 6) / seq(1, 6)  # 48/1, 48/2 ... 48/12
loadfile <- "Robjs/liver_kidney_atger_nestle/dat.complex.all_T.Robj"
if (file.exists(loadfile)){
  load(loadfile)
} else {
  
  library(parallel)
  dat.complexes <- mclapply(periods, function(period, dat.long){
    omega <- 2 * pi / period
    dat.tmp <- dat.long %>%
      group_by(gene, tissue) %>%
      do(ProjectToFrequency2(., omega, add.tissue=TRUE))
    dat.tmp$period <- period
    return(dat.tmp)
  }, dat.long = dat.long, mc.cores = length(periods))
  
  
  dat.complex.all_T <- do.call(rbind, dat.complexes)
  outfcomp <- "Robjs/liver_kidney_atger_nestle/dat.complex.all_T.bugfixed.Robj"
  if (!file.exists(outfcomp)) save(dat.complex.all_T, file = outfcomp)
  rm(dat.complexes)
}
outffreq <- "Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj"
if (!file.exists(outffreq)) save(dat.freq, file = outffreq)




jmeth <- "BIC"
jmeth <- "g=1001"
i <- 1

fits.count <- subset(fits.long.filt, method == jmeth & model != "") %>% group_by(model) %>% summarise(model.count = length(model))
fits.count <- fits.count[order(fits.count$model.count, decreasing = TRUE), ]
fits.count <- subset(fits.count, model.count > 190)  # no underdetermination
jmodel.lst <- as.character(fits.count$model)

N.f <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/promoters_filtered_by_gene_liver_kidney_novel_modules/sitecount_matrix_geneids.Liver_SV129.g=1001.mat"

N <- read.table(N.f, header=TRUE)
jmotif.order <- "RORA.p2"
N.ordered <- N[order(N[[jmotif.order]], decreasing=TRUE), ]
