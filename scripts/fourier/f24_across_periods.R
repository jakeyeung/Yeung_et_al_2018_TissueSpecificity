# f24_across_periods.R
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
# # across all genes
# start <- Sys.time()
# dat.fit.periods.genome_wide <- FitRhythmicScanPeriods(dat.long, periods, cores = 20)
# print(Sys.time() - start)
# head(dat.fit.periods.genome_wide)
# save(dat.fit.periods.genome_wide, file = "Robjs/dat.fit.scan_periods.genome_wide.Robj")
# 2.5 hours later...

# load(file = "Robjs/dat.fit.scan_periods.Robj")
# dat.fit.periods <- subset(dat.fit.periods, gene %in% clockgenes)
# 
# dat.fit.periods$chi.sqr <- dat.fit.periods$ssq.residuals / dat.fit.periods$variance
# 
# ggplot(dat.fit.periods, aes(x = period, y = ssq.residuals, colour = gene, group = gene)) + geom_line() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")
# ggplot(subset(dat.fit.periods, ! gene %in% "Asb12"), aes(x = period, y = chi.sqr, colour = gene, group = gene)) + geom_line() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")
# 
# 
# # Get minimum per gene ----------------------------------------------------
# 
# dat.fit.periods.min <- dat.fit.periods %>%
#   group_by(gene, tissue) %>%
#   do(GetMinPeriodSsqResiduals(.))
# 
# ggplot(dat.fit.periods.min, aes(x = period, y = chi.sqr, colour = gene)) + geom_point() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")
# 
# ggplot(dat.fit.periods.min, aes(x = period)) + geom_histogram() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")
# ggplot(dat.fit.periods.min, aes(x = period)) + geom_histogram() + geom_vline(xintercept=24, linetype="dotted")

# dat.fit.periods.mean <- dat.fit.periods %>%
#   group_by(tissue, period) %>%
#   summarise(mean.ssq.residuals = mean(ssq.residuals / variance), sd.ssq.residuals = sd(ssq.residuals / variance))

# ggplot(subset(dat.fit.periods.mean, tissue %in% c("BS", "Cere", "Hypo")), aes(x = period, y = mean.ssq.residuals)) + geom_line() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")
# ggplot(dat.fit.periods.mean, aes(x = period, y = mean.ssq.residuals)) + geom_line() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")

# dat.fit.periods.min <- dat.fit.periods.mean %>%
#   group_by(tissue) %>%
#   do(GetMinPeriod(.))
# dat.fit.periods.min

# ggplot(dat.fit.periods.min, aes(x = tissue, y = period)) + geom_bar(stat = "identity")
# ggplot(dat.fit.periods.min, aes(x = tissue, y = period)) + geom_bar(stat = "identity")


# Bayesian linear regression ----------------------------------------------

# y <- subset(dat.long, gene == "Nr1d1" & tissue == "Liver" & experiment == "array")$exprs
# x <- subset(dat.long, gene == "Nr1d1" & tissue == "Liver" & experiment == "array")$time
# plot(x, y)
