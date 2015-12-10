# 2015-12-10
# Jake Yeung
# same as dhs_peak_gene_body.R but across all motifs.

library(hash)
library(dplyr)
library(ggplot2)
library(reshape2)

source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
# Load objs -------------------------------------------------------------

load("Robjs/N.long.liver_genes.all_motifs.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.long.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)


# Functions ---------------------------------------------------------------

TakeIfSame <- function(x){
  if (length(unique(x)) > 1){
    return(NA)
  } else {
    return(mean(x))
  }
}

# Get cutoffs ---------------------------------------------------------

liver.genes <- as.character(subset(fits.best, model == "Liver")$gene)
flat.genes <- as.character(subset(fits.best, model == "")$gene)

# cutoff for the DHS

# shift by 1 log2 unit the fit looks better
log.shift <- 1
S.tissuecutoff$cutoff.adj <- 2^(log2(S.tissuecutoff$cutoff) + log.shift)

# cool
ggplot(S.long[sample(x = 1:nrow(S.long), size = 0.01 * nrow(S.long)), ], aes(x = log2(signal + pseudo))) + geom_density() + 
  geom_vline(aes(xintercept = log2(cutoff.adj)), data = S.tissuecutoff) + facet_wrap(~tissue)

cutoffs.tiss <- hash(as.character(S.tissuecutoff$tissue), as.numeric(S.tissuecutoff$cutoff.adj))


# Get peaks by cutoffs ----------------------------------------------------


S.sub <- subset(S.long, gene %in% liver.genes)


S.sub$signal.cut <- mapply(function(s, tiss){
  cutoff.tiss <- cutoffs.tiss[[tiss]]
  if (s >= cutoff.tiss){
    return(1)
  } else {
    return(0)
  }
}, S.sub$signal, as.character(S.sub$tissue))


start <- Sys.time()
dist.filt <- 1000
S.collapse <- subset(S.sub, dist <= dist.filt) %>%
  group_by(gene, peak) %>%
  do(CollapseDat(., indx = 4, tissue = "Liver", non.tissue = "Flat", flat.style = "normal"))
print(Sys.time() - start)


# Do enrichment -----------------------------------------------------------

N.sub.all <- subset(N.long, peak %in% S.collapse$peak)

RunFisher <- function(dat, S, sitecounts.hash){
  sitecounts.hash <- hash(as.character(dat$peak), dat$sitecount)
  S$sitecount <- sapply(S$peak, AssignSitecount, sitecounts.hash)
  test <- FisherTestSitecounts(dat = S, cutoff = 0.5, sitecount.col = "sitecount", 
                       model.col = "peak.type", show.table=FALSE)
  return(test)
}

tests.motif <- N.sub.all %>%
  group_by(motif) %>%
  do(RunFisher(., S.collapse, sitecounts.hash))

# plot volcano
ggplot(tests.motif, aes(x = odds.ratio, y = -log10(p.value), label = motif)) + geom_text()

# # single motif
# jmotif <- "ONECUT1,2.p2"
# N.sub <- subset(N.sub.all, motif == jmotif)
# dim(N.sub)
# sitecounts.hash <- hash(as.character(N.sub$peak), N.sub$sitecount)
# S.collapse$sitecount <- sapply(S.collapse$peak, AssignSitecount, sitecounts.hash)
# FisherTestSitecounts(dat = S.collapse, cutoff = 0.5, sitecount.col = "sitecount", model.col = "peak.type", show.table=TRUE)


# Do between liver and flat genes -----------------------------------------

S.sub.lf <- subset(S.long, gene %in% c(liver.genes, flat.genes))

S.sub.lf$signal.cut <- mapply(function(s, tiss){
  cutoff.tiss <- cutoffs.tiss[[tiss]]
  if (s >= cutoff.tiss){
    return(1)
  } else {
    return(0)
  }
}, S.sub.lf$signal, as.character(S.sub.lf$tissue))


# remove peaks assigned to multiple genes
S.sub.lf <- subset(S.sub.lf, dist <= 10000) %>%
  group_by(peak, tissue) %>%
  # summarise(signal.cut = TakeIfSame(signal.cut))
  summarise(signal.cut = mean(signal.cut))

S.sub.lf.mat <- dcast(data = subset(S.sub.lf), formula = peak ~ tissue, value.var = "signal.cut")
rownames(S.sub.lf.mat) <- S.sub.lf.mat$peak; S.sub.lf.mat$peak <- NULL

# find liver-specific peaks
S.is.livpeaks <- apply(S.sub.lf.mat, 1, function(signal.cut){
  liver.indx <- 4
  liver.sig <- signal.cut[liver.indx]
  other.sig <- signal.cut[-liver.indx]
  # normal
  if (liver.sig == 1 & max(other.sig) == 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
})

S.sub.lf.mat.filt <- S.sub.lf.mat[S.is.livpeaks, ]

# annotate peaks to either Liver or Flat model
S.sub.l <- as.character(subset(S.long, gene %in% liver.genes)$peak)
S.sub.f <- as.character(subset(S.long, gene %in% flat.genes)$peak)

jkey <- c(S.sub.l, S.sub.f)
jval <-  c(rep("Liver", length(S.sub.l)), rep("Flat", length(S.sub.f)))
peak.hash <- hash(jkey, jval)

# annotate
N.sub <- subset(N.RORA, gene %in% c(liver.genes, flat.genes))

N.sub$model <- sapply(as.character(N.sub$peak), function(p) peak.hash[[p]])

# test one motif
FisherTestSitecounts(subset(N.sub, motif == "RORA.p2"), cutoff = 0.5, sitecount.col = "sitecount", model.col = "model", show.table = TRUE)

# Run across motifs
tests.motif <- N.sub %>%
  group_by(motif) %>%
  do(RunFisher(., S.collapse, sitecounts.hash))

# crashes
# start <- Sys.time()
# dist.filt <- 1000
# S.collapse.lf <- subset(S.sub.lf, dist <= dist.filt) %>%
#   group_by(gene, peak) %>%
#   do(CollapseDat(., indx = 4, tissue = "Liver", non.tissue = "Flat", flat.style = "normal"))
# print(Sys.time() - start)
