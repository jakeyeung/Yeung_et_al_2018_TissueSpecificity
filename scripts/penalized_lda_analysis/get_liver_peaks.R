# 2016-02-17
# get_liver_peaks.R
# be super stringent in identifying liver peaks

library(ggplot2)
library(dplyr)
library(hash)

setwd("/home/yeung/projects/tissue-specificity")

load("Robjs/S.long.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)


# Functions ---------------------------------------------------------------

source("scripts/functions/SitecountsFunctions.R")

# Find cutoff -------------------------------------------------------------

# shift by 1 log2 unit the fit looks better
log.shift <- 2.5
S.tissuecutoff$cutoff.adj <- 2^(log2(S.tissuecutoff$cutoff) + log.shift)

cutoff.adj <- S.tissuecutoff$cutoff.adj
cutoff.lower <- S.tissuecutoff$cutoff

# show upper and lower limit
pseudo <- 1e-3
ggplot(S.long[sample(x = 1:nrow(S.long), size = 0.01 * nrow(S.long)), ], aes(x = log2(signal + pseudo))) + geom_density() + 
  facet_wrap(~tissue) +
  geom_vline(aes(xintercept = log2(cutoff.adj)), data = S.tissuecutoff, colour = "blue") + 
  geom_vline(aes(xintercept = log2(cutoff)), data = S.tissuecutoff, colour = "red")

cutoffs.tiss.upper <- hash(as.character(S.tissuecutoff$tissue), as.numeric(S.tissuecutoff$cutoff.adj))
cutoffs.tiss.lower <- hash(as.character(S.tissuecutoff$tissue), as.numeric(S.tissuecutoff$cutoff))


# Get liver-specific peaks  -----------------------------------------------

livflat.genes <- subset(fits.best, model %in% c("Liver", ""))$gene

jstart <- Sys.time()  # 4 min 
S.sub.collapse <- subset(S.long, gene %in% livflat.genes)
S.sub.collapse$is.lower <- mapply(function(tiss, sig, cutoffs.lower){
  IsSignalLower(tiss, sig, cutoffs.lower)
  }, S.sub.collapse$tissue, S.sub.collapse$signal, MoreArgs = list(cutoffs.lower = cutoffs.tiss.lower))
S.sub.collapse$is.upper <- mapply(function(tiss, sig, cutoffs.upper){
  IsSignalUpper(tiss, sig, cutoffs.upper)
  }, S.sub.collapse$tissue, S.sub.collapse$signal, MoreArgs = list(cutoffs.upper = cutoffs.tiss.upper))
print(Sys.time() - jstart)

# filter out peaks that are above threshold in liver but below threshold in other tissues

# jtmp <- subset(fits.best, model == "Liver" & amp.avg > 0.2)
# top.genes <- as.character(jtmp[order(jtmp$amp.avg, decreasing=TRUE), ]$gene)
jtmp <- subset(fits.best, model == "")
# too many genes it crashes take top 500
frac <- 0.7
top.n <- round(nrow(jtmp) * 0.7)
top.genes <- as.character(jtmp[order(jtmp$amp.avg, decreasing=TRUE), ]$gene[1:top.n])
peaks <- unique(as.character(subset(S.sub.collapse, gene %in% top.genes)$peak))

# crashes, do subset of S.sub.collapse instead
S.sub.collapse.filt <- subset(S.sub.collapse, peak %in% peaks) %>%
# S.sub.collapse.filt <- S.sub.collapse %>%
  group_by(peak, gene) %>%
  do(IsTissueSpecificLong(.))

print(subset(S.sub.collapse.filt, is.tiss.spec == TRUE))
print(length(unique(subset(S.sub.collapse.filt, is.tiss.spec == TRUE)$gene)))

save(S.sub.collapse.filt, file = "Robjs/S.flat.liverspec.peaks.top70percent.Robj")
