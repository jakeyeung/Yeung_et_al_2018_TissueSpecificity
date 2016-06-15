# 2016-06-12
# Jake Yeung 

rm(list=ls())
library(dplyr)
library(ggplot2)
setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FitRhythmicAcrossPeriods.R")
source("scripts/functions/LiverKidneyFunctions.R")

# Load --------------------------------------------------------------------

dat <- LoadLivKid()
foutpath <- "Robjs/liver_kidney/dat.fit.periods.genome_wide.Robj"
load(foutpath, verbose = T)

dat.fit.periods.genome_wide.min <- dat.fit.periods.genome_wide %>%
  group_by(gene, tissue) %>%
  filter(period == period[which.min(ssq.residuals)])

xscale_periods <- seq(9, 31, 2)
plot.periods.all <- ggplot(dat.fit.periods.genome_wide.min, aes(x = period)) + 
  geom_histogram(binwidth = 0.5) + 
  geom_vline(xintercept=24, linetype="dotted") + 
  scale_x_continuous(breaks=xscale_periods) + 
  xlab("Period with best fit [h]") + ylab("Number of genes") +
  theme_bw(24) + theme(aspect.ratio=1,
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       axis.text.x=element_text(size=11))
print(plot.periods.all)
print(plot.periods.all + facet_wrap(~tissue))


