# f24_across_periods.R
# Jake Yeung
# 2015-09-22
# library("devtools")
# dev_mode()
# 
# install("~/projects/f24")  # use jake branch
# library(f24.R2.cycling)

library(dplyr)
setwd("/home/yeung/projects/tissue-specificity")

FitRhythmicScanPeriods <- function(dat.long, periods, cores = 51){
  # likely a filtered dat.long set
  library(parallel)
  #   periods <- list(periods)
  dat.fitrhyth.period <- mclapply(periods, function(p){
    dat.fit <- dat.long %>%
      group_by(gene, tissue) %>%
      do(FitRhythmic(., T = p, get.residuals=TRUE))
    dat.fit$period = p
    return(dat.fit)
  }, mc.cores = 64)
  dat.fitrhyth <- do.call(rbind, dat.fitrhyth.period)
  return(dat.fitrhyth)
}

GetMinPeriod <- function(dat){
  return(subset(dat, mean.ssq.residuals == min(dat$mean.ssq.residuals)))
}



# Source ------------------------------------------------------------------

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")

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
dat.fit.periods <- FitRhythmicScanPeriods(subset(dat.long, gene %in% clockgenes), periods, cores = 51)
head(dat.fit.periods)
save(dat.fit.periods, file = "Robjs/dat.fit.scan_periods.Robj")

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
