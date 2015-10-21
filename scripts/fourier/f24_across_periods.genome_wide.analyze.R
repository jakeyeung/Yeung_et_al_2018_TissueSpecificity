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


# Perform f24 -------------------------------------------------------------

# on a single gene
jgene <- "Nr1d1"
jtissue <- "Liver"
clockgenes <- GetClockGenes()

period <- 24
dat.fit <- subset(dat.long, gene %in% clockgenes) %>%
  group_by(gene, tissue) %>%
  do(FitRhythmic(., T = period, get.residuals = TRUE))
dat.fit
# dat.fit <- FitRhythmic(dat.long)

p.min <- 20
p.max <- 30
periods <- seq(p.min, p.max, by = 0.1)
# across all genes
# start <- Sys.time()
# dat.fit.periods.genome_wide <- FitRhythmicScanPeriods(dat.long, periods, cores = 20)
# print(Sys.time() - start)
# head(dat.fit.periods.genome_wide)
# save(dat.fit.periods.genome_wide, file = "Robjs/dat.fit.scan_periods.genome_wide.Robj")
# 2.5 hours later...

load("Robjs/dat.fit.scan_periods.genome_wide.Robj", verbose=T)

clockgenes <- GetClockGenes()

dat.fit.periods.genome_wide$chi.sqr <- dat.fit.periods.genome_wide$ssq.residuals / dat.fit.periods.genome_wide$variance
dat.fit.periods.genome_wide.min <- dat.fit.periods.genome_wide %>%
  group_by(gene, tissue) %>%
  filter(period == period[which.min(ssq.residuals)])

# dat.fit.periods.min.check <- subset(dat.fit.periods.genome_wide.min, gene %in% unique(dat.fit.periods.min$gene))  # from .min from clockgenes
dat.fit.periods.sub <- subset(dat.fit.periods.genome_wide.min, amp > 1 & pval < 1e-5)
  
ggplot(subset(dat.fit.periods.sub, tissue != "WFAT"), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.genome_wide.min$period))/100) + facet_wrap(~tissue)
ggplot(subset(dat.fit.periods.sub, tissue != "WFAT"), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.genome_wide.min$period))/100) + geom_vline(xintercept=24, linetype="dotted")


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
