FitRhythmicScanPeriods <- function(dat.long, periods, cores = 51){
  # likely a filtered dat.long set
  source("~/projects/tissue-specificity/scripts/functions/FitRhythmic.R")
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

GetMinPeriodSsqResiduals <- function(dat){
  return(subset(dat, chi.sqr == min(dat$chi.sqr)))
}