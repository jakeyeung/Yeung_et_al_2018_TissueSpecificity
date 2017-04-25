# 2017-04-20
# Jake Yeung
# Use Weight cutoff and then run three way contingency table
# What makes FOX CUX ONECUT interesting? Or are others more interesting?

rm(list=ls())

library(dplyr)
library(ggplot2)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

RunPoissonModel.test <- function(dat){
  mod1 <- glm(freq ~ model + motif1 * motif2, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.JI <- function(dat){
  mod1 <- glm(freq ~ model + motif1 + motif2 + motif1 * motif2, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.sing <- function(dat){
  mod1 <- glm(freq ~ model + motif1, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.cnd <- function(dat){
  mod1 <- glm(freq ~ model + motif1 + motif2 + model * motif1 + model * motif2, data=dat, family=poisson())
  return(mod1)
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}

RunPoissonModel.cnd2 <- function(dat){
  mod1 <- glm(freq ~ (motif1 + motif2) * model, data=dat, family=poisson())
  chisq.pval <- 1 - pchisq(mod1$deviance, mod1$df.residual)
  return(data.frame(deviance = mod1$deviance, pval = chisq.pval))
}


# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)

K <- 200
weight.cutoff <- 0.8
# inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.withNmatallNmatfreqs.RemoveZeroCounts.Robj"
inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.", K, ".weight.", weight.cutoff, ".withNmatallNmatfreqs.RemoveZeroCounts.Robj")

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.200.weight.0.8.MergePeaks.TRUE.withNmatallNmatfreqs.RemoveZeroCounts.Robj"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts.LiverPeaksvsBoth.Robj"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts.LiverPeaksvsBoth.Robj"
inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts2.Robj"

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.nullmodel.CI.withNmatallNmatfreqs.RemoveZeroCounts.Robj"


inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.nullmodel.CI.withNmatallNmatfreqs.RemoveZeroCounts.Robj"


inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts2.Robj"

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.3.K.300.weight.0.8.MergePeaks.FALSE.nullmodel.CI.withNmatallNmatfreqs.RemoveZeroCounts.Robj"

inf <- "/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.3.K.300.weight.0.8.MergePeaks.FALSE.withNmatallNmatfreqs.RemoveZeroCounts2.Robj"


K <- 300
weight <- 0.8
nb.models <- 3
inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.", nb.models, ".K.", K, ".weight.", weight, ".MergePeaks.FALSE.nullmodel.JI.withNmatallNmatfreqs.RemoveZeroCounts.Robj")

load(inf, v=T)  # load also Nmat and Nmatfreqs
N.mat.freqs.all <- N.mat.freqs

if (length(unique(N.mat.freqs.all$model)) == 3){
  N.mat.freqs <- subset(N.mat.freqs, model %in% c("flat", "rhyth"))
}


# Do fit on pairs: but merge flat and kidney  -----------------------------

N.mat.freqs.all$model.merge <- sapply(N.mat.freqs.all$model, function(m) ifelse(m == "rhyth", "rhyth", "flat"))

N.mat.freqs.merged <- N.mat.freqs.all %>%
  group_by(model.merge, pair, motif1, motif2) %>%
  summarise(freq = sum(freq))

N.mat.freqs.merged <- dplyr::rename(N.mat.freqs.merged, "model" = model.merge)

fits.merged <- subset(N.mat.freqs.merged) %>%
  group_by(pair) %>%
  do(RunPoissonModel.JI(.))

fits <- fits.merged

save(fits, N.mat.all, N.mat.freqs, file = paste0("Robjs/three_way_cooccurence/three.way.cooccurrence.LiverVsKidFLATMerged.Robj"))

