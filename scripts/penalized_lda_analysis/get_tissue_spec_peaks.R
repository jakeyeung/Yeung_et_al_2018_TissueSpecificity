# 2016-02-19
# get_tissue_spec_peaks.R
# Find tissue specific peaks within liver-specific rhythmic genes

library(hash)

# Load --------------------------------------------------------------------

load("Robjs/S.long.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)

# Get liver genes ---------------------------------------------------------

ampmin <- 0.25
liv.genes <- unique(as.character(subset(fits.best, model == "Liver" & amp.avg > ampmin)$gene))

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

jstart <- Sys.time()  # 4 min 
S.sub.collapse <- subset(S.long, gene %in% liv.genes)

S.sub.collapse$is.upper <- apply(S.sub.collapse, 1, function(row){
  tiss <- row[5]
  sig <- row[6]
  IsSignalUpper(tiss, sig, cutoffs.upper)
})

S.sub.collapse$is.upper <- mapply(function(tiss, sig, cutoffs.upper){
  IsSignalUpper(tiss, sig, cutoffs.upper)
}, S.sub.collapse$tissue, S.sub.collapse$signal, MoreArgs = list(cutoffs.upper = cutoffs.tiss.upper))
print(Sys.time() - jstart)
