library(dplyr)
library(parallel)
library(ggplot2)

rm(list=ls())

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/fits.liver_kidney.Robj", v=T)



source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/NcondsFunctions.R")

# get best model from each n.params
fits.best.livkid <- fits.all.long %>%
  group_by(gene) %>%
  filter(weight.raw == min(weight.raw))

fits.best.livkid$model <- mapply(FilterModelByAmp, fits.best.livkid$model, fits.best.livkid$param.list, MoreArgs = list(amp.cutoff = 0.15))



# phase distribution
ggplot(subset(fits.best.livkid, model == "Liver;Kidney"), aes(x = phase.maxdiff)) + geom_histogram(bins = 40)
ggplot(subset(fits.best.livkid, model == "Liver;Kidney"), aes(x = amp.avg)) + geom_histogram(bins = 40)
ggplot(subset(fits.best.livkid, model == "Liver"), aes(x = amp.avg)) + geom_histogram(bins = 40)
ggplot(subset(fits.best.livkid, model == "Kidney"), aes(x = amp.avg)) + geom_histogram(bins = 40)
ggplot(subset(fits.best.livkid, model == "Liver,Kidney"), aes(x = amp.avg)) + geom_histogram(bins = 40)
