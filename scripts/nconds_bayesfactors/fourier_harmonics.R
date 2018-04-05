# 2016-06-12
# Jake Yeung
# fourier_harmonics.R

rm(list=ls())

source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/SvdFunctions.R")

# Load --------------------------------------------------------------------

dat <- LoadLivKid()

periods <- rep(48, 12) / seq(1, 12)  # 48/1, 48/2 ... 48/12

library(parallel)
dat.complexes <- lapply(periods, function(period, dat, add.entropy.method){
  dat.comp <- TemporalToFrequencyDatLong(dat = dat, period = period, add.entropy.method = add.entropy.method)
  dat.comp$period <- period
  return(dat.comp)
  # }, dat = dat.long, add.entropy.method = FALSE)
}, dat = dat, add.entropy.method = FALSE)
dat.complex.all_T <- do.call(rbind, dat.complexes)


# Variance of each s  -----------------------------------------------------

dat.var.s <- dat.complex.all_T %>%
  group_by(period, tissue) %>%
  summarise(sum_sqr_mod = sum(Mod(exprs.transformed) ^ 2)) %>%
  mutate(period.factor = signif(period, digits = 3))

dat.var.s$period.factor <- factor(dat.var.s$period.factor, 
                                  levels = sort(unique(dat.var.s$period.factor), decreasing = TRUE))
dat.var.s <- dat.var.s %>%
  group_by(tissue) %>%
  mutate(sum_sqr_mod.norm = sum_sqr_mod / sum(sum_sqr_mod))

ggplot(dat.var.s, aes(x = period.factor, y = sum_sqr_mod)) +  geom_bar(stat = "identity") + facet_wrap(~tissue) + 
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) + xlab("Period [h]") + ylab("Amplitude squared")
ggplot(dat.var.s, aes(x = period.factor, y = sum_sqr_mod.norm)) +  geom_bar(stat = "identity") + facet_wrap(~tissue) + 
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) + xlab("Period [h]") + ylab("Normalized amplitude squared")
 

