# 2016-06-10
# Jake Yeung
# downstream_analysis.R

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")

# Load --------------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/Robjs/bayes_factors"
inf <- list.files(indir)

fits.long <- expandingList()
for (f in inf){
  method <- strsplit(f, split = "\\.")[[1]][[5]]
  fpath <- file.path("Robjs/bayes_factors", f)
  load(fpath)
  fits.all.long$method <- method
  fits.long$add(fits.all.long)
}
fits.long <- fits.long$as.list()
fits.long <- do.call(rbind, fits.long)


# Summarize ---------------------------------------------------------------

fits.long.filt <- fits.long %>%
  group_by(gene, method) %>%
  filter(weight.raw == min(weight.raw))

dat <- LoadLivKid()

PlotGeneAcrossTissues(subset(dat, gene == "Cry1"))


# Add amplitude and phase info --------------------------------------------

fits.long.filt$amp <- sapply(fits.long.filt$param.list, function(p) GetAvgAmpFromParams(params = p, by.model = FALSE))
fits.long.filt$phase <- sapply(fits.long.filt$param.list, function(p) GetPhaseFromParams(params = p, by.model = FALSE))
fits.long.filt$phase.diff <- sapply(fits.long.filt$param.list, function(p) GetMaxPhaseDiffFromParams(params = p, by.model=FALSE))

# Count -------------------------------------------------------------------

fits.sum <- fits.long.filt %>%
  group_by(method, model) %>%
  summarise(count = length(model))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(fits.sum, aes(x = model, fill = method, y = count)) + geom_bar(stat = "identity", position = "dodge", width=0.5) + theme_bw(24) + scale_fill_manual(values=cbPalette)

# Core-clock genes --------------------------------------------------------

jsub <- subset(fits.long.filt, method == "zf")
print(jsub[order(jsub$weight, decreasing = TRUE), ])


# Show kidney-only genes --------------------------------------------------

jsub <- subset(fits.long.filt, method == "zf" & model == "Kidney")
print(jsub[order(jsub$weight, decreasing = TRUE), ])
