# f24_across_periods.R
# CAN RUN THIS ON SERVER AS IS TAKES 2.5 HOURS
# Jake Yeung
# 2015-09-22
# library("devtools")
# dev_mode()
# 
# install("~/projects/f24")  # use jake branch
# library(f24.R2.cycling)

library(dplyr)
library(ggplot2)
setwd("/home/yeung/projects/tissue-specificity")



# Source ------------------------------------------------------------------

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FitRhythmicAcrossPeriods.R")

# Load data ---------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
# load("Robjs/dat.long.Robj")


# load("Robjs/dat.fit.scan_periods.genome_wide.Robj", verbose=T)
load("/home/yeung/projects/tissue-specificity/Robjs/dat.fit.scan_periods.genome_wide.10_to_30.Robj", verbose=T)
dat.fit.periods.genome_wide.10_to_30 <- dat.fit.periods.genome_wide
load("/home/yeung/projects/tissue-specificity/Robjs/dat.fit.scan_periods.genome_wide.5_10.Robj", verbose=T)
dat.fit.periods.genome_wide <- rbind(dat.fit.periods.genome_wide, dat.fit.periods.genome_wide.10_to_30)

dat.fit.periods.genome_wide.min <- dat.fit.periods.genome_wide %>%
  group_by(gene, tissue) %>%
  filter(period == period[which.min(ssq.residuals)])

# dat.fit.periods.min.check <- subset(dat.fit.periods.genome_wide.min, gene %in% unique(dat.fit.periods.min$gene))  # from .min from clockgenes
dat.fit.periods.sub <- subset(dat.fit.periods.genome_wide.min, amp > 0.5 & pval < 1e-5)
  
ggplot(subset(dat.fit.periods.sub, tissue != "WFAT"), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.genome_wide.min$period))/100) + facet_wrap(~tissue)
ggplot(subset(dat.fit.periods.sub, tissue != "WFAT"), aes(x = period)) + 
  geom_histogram(binwidth = diff(range(dat.fit.periods.genome_wide.min$period))/100) + 
  geom_vline(xintercept=24, linetype="dotted") + geom_vline(xintercept=12, linetype="dotted")

tiss <- "BFAT"; gen <- "Myf6"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Adr"; gen <- "Elovl3"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "BFAT"; gen <- "Ampd1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "BFAT"; gen <- "Ampd1"; exper="rnaseq"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Liver"; gen <- "Thrsp"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Mus"; gen <- "Tnni1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)
