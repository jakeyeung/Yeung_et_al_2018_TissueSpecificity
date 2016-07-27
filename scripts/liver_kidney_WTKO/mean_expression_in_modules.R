# 2016-07-27
# Jake Yeung
# Are modules differentially expressed? Liver and Kidney

rm(list=ls())

source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotFunctions.R")


library(dplyr)
library(ggplot2)

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)

dat.orig <- dat.long
dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)


# Get mean ----------------------------------------------------------------

dat.mean <- dat.long %>%
  group_by(gene, tissue, geno) %>%
  summarise(exprs.mean = mean(exprs))

jmeth <- "g=1001"
fits.count <- subset(fits.long.filt, method == jmeth & model != "") %>%
  group_by(model) %>%
  summarise(count = length(gene)) %>%
  arrange(desc(count))

jmods <- as.character(subset(fits.count, count > 190)$model)

jmod <- "Kidney_SV129"
for (jmod in jmods){
  print(jmod)
  genes <- as.character(subset(fits.long.filt, model == jmod & method == jmeth)$gene)
  print(PlotMeanExprsOfModel(dat.mean, genes = genes, jmodel = jmod, sorted = TRUE, avg.method = "mean"))
}
