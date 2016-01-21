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

N.long <- subset(N.long, dist <= 1000)
S.long <- subset(S.long, dist <= 1000)

# Functions ---------------------------------------------------------------

FilterByRange <- function(N, ucsc_coord, cutoff = 0.5){
  # parse chromo, start, end from copy/paste of ucsc_coord
  jchromo <- strsplit(ucsc_coord, ":")[[1]][[1]]
  startend <- strsplit(ucsc_coord, ":")[[1]][[2]]
  jstart <- strsplit(startend, "-")[[1]][[1]]
  jend <- strsplit(startend, "-")[[1]][[2]]
  # parse out commas and turn into int
  jstart <- as.numeric(gsub(",", "", x = jstart))
  jend <- as.numeric(gsub(",", "", x = jend))
  sub <- subset(N, chromo == jchromo & start > jstart & end < jend & sitecount > cutoff)
  return(sub[order(sub$sitecount, decreasing = TRUE), ])
}



# Get cutoffs ---------------------------------------------------------

liver.genes <- as.character(subset(fits.best, model == "Liver")$gene)
flat.genes <- as.character(subset(fits.best, model == "")$gene)

# cutoff for the DHS

# shift by 1 log2 unit the fit looks better
log.shift <- 1.5
S.tissuecutoff$cutoff.adj <- 2^(log2(S.tissuecutoff$cutoff) + log.shift)

# cool
pseudo <- 1e-3
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

# remove peaks assigned to multiple genes
# test.out <- test %>%
#   group_by(peak, tissue) %>%
#   summarise(signal.cut = unique(signal.cut),
#             gene = paste(gene, collapse = ";"),
#             dist = unique(dist),
#             signal.cut = unique(signal.cut))

S.sub <- S.sub %>%
  group_by(peak, tissue) %>%
  summarise(signal.cut = unique(signal.cut),
            gene = paste(gene, collapse = ";"),
            dist = unique(dist),
            signal.cut = unique(signal.cut))

start <- Sys.time()
dist.filt <- 1000
S.collapse <- subset(S.sub, dist <= dist.filt) %>%
  group_by(gene, peak) %>%
  do(CollapseDat(., indx = 4, tissue = "Liver", non.tissue = "Flat", flat.style = "normal"))
print(Sys.time() - start)


# Do enrichment -----------------------------------------------------------

N.sub.all <- subset(N.long, peak %in% S.collapse$peak)


tests.motif <- N.sub.all %>%
  group_by(motif) %>%
  do(RunFisherDHS(., S.collapse))

# plot volcano
ggplot(tests.motif, aes(x = odds.ratio, y = -log10(p.value), label = motif)) + geom_text()

# # single motif
# jmotif <- "ONECUT1,2.p2"
# N.sub <- subset(N.sub.all, motif == jmotif)
# dim(N.sub)
# sitecounts.hash <- hash(as.character(N.sub$peak), N.sub$sitecount)
# S.collapse$sitecount <- sapply(S.collapse$peak, AssignSitecount)
# FisherTestSitecounts(dat = S.collapse, cutoff = 0.5, sitecount.col = "sitecount", model.col = "peak.type", show.table=TRUE)


# Add flat genes ----------------------------------------------------------

S.sub.flat <- subset(S.long, gene %in% flat.genes)

S.sub.flat$signal.cut <- mapply(function(s, tiss){
  cutoff.tiss <- cutoffs.tiss[[tiss]]
  if (s >= cutoff.tiss){
    return(1)
  } else {
    return(0)
  }
}, S.sub.flat$signal, as.character(S.sub.flat$tissue))

S.sub.flat <- S.sub.flat %>%
  group_by(peak, tissue) %>%
  summarise(signal.cut = unique(signal.cut),
            gene = paste(gene, collapse = ";"),
            dist = unique(dist),
            signal.cut = unique(signal.cut))

# takes ~90 minutes
start <- Sys.time()
dist.filt <- 1000
S.collapse.flat <- subset(S.sub.flat, dist <= dist.filt) %>%
  group_by(gene, peak) %>%
  do(CollapseDat(., indx = 4, tissue = "Liver", non.tissue = "Flat", flat.style = "normal"))
print(Sys.time() - start)

# label different models
S.collapse.flat$model <- "Flat"
S.collapse$model <- "Liver"

# merge together
S.collapse <- rbind(S.collapse, S.collapse.flat)

# filter for liver peaks
S.collapse <- subset(S.collapse, peak.type == "Liver")

# load up sitecounts for liver and flat genes
load("Robjs/N.long.flatliver_genes.all_motifs.dist1000.Robj", verbose=T)

N.long.all <- subset(N.long.all, peak %in% S.collapse$peak)

tests.motif.livflat <- N.long.all  %>%
  group_by(motif) %>%
  do(RunFisherDHS(., S.collapse, jmodel.col = "model"))

# plot volcano

RunFisherDHS(subset(N.long.all, motif == "HIF1A.p2"), S.collapse, jmodel.col = "model", jshow.table = TRUE)

ggplot(tests.motif.livflat, aes(x = odds.ratio, y = -log10(p.value), label = motif)) + geom_text()


# Test on Celsr1 ----------------------------------------------------------

load("Robjs/S.collapse.liver.Robj", verbose=T)

S.test <- subset(S.collapse, gene == "Celsr1")
subset(S.test, peak.type == "Liver")

jpeak <- "chr15:85961695-85962195"
jpeak <- "chr15:85962695-85963195"
subset(N.long.all, peak == jpeak & sitecount > 0.5)


# # Do between liver and flat genes -----------------------------------------
# 
# S.sub.lf <- subset(S.long, gene %in% c(liver.genes, flat.genes))
# 
# S.sub.lf$signal.cut <- mapply(function(s, tiss){
#   cutoff.tiss <- cutoffs.tiss[[tiss]]
#   if (s >= cutoff.tiss){
#     return(1)
#   } else {
#     return(0)
#   }
# }, S.sub.lf$signal, as.character(S.sub.lf$tissue))
# 
# 

# 
# S.sub.lf.mat <- dcast(data = subset(S.sub.lf), formula = peak ~ tissue, value.var = "signal.cut")
# rownames(S.sub.lf.mat) <- S.sub.lf.mat$peak; S.sub.lf.mat$peak <- NULL
# 
# # find liver-specific peaks
# S.is.livpeaks <- apply(S.sub.lf.mat, 1, function(signal.cut){
#   liver.indx <- 4
#   liver.sig <- signal.cut[liver.indx]
#   other.sig <- signal.cut[-liver.indx]
#   # normal
#   if (liver.sig == 1 & max(other.sig) == 0){
#     return(TRUE)
#   } else {
#     return(FALSE)
#   }
# })
# 
# S.sub.lf.mat.filt <- S.sub.lf.mat[S.is.livpeaks, ]
# 
# # annotate peaks to either Liver or Flat model
# S.sub.l <- as.character(subset(S.long, gene %in% liver.genes)$peak)
# S.sub.f <- as.character(subset(S.long, gene %in% flat.genes)$peak)
# 
# jkey <- c(S.sub.l, S.sub.f)
# jval <-  c(rep("Liver", length(S.sub.l)), rep("Flat", length(S.sub.f)))
# peak.hash <- hash(jkey, jval)
# 
# # annotate
# N.sub <- subset(N.RORA, gene %in% c(liver.genes, flat.genes))
# 
# N.sub$model <- sapply(as.character(N.sub$peak), function(p) peak.hash[[p]])
# 
# # test one motif
# FisherTestSitecounts(subset(N.sub, motif == "RORA.p2"), cutoff = 0.5, sitecount.col = "sitecount", model.col = "model", show.table = TRUE)
# 
# # Run across motifs
# tests.motif <- N.sub %>%
#   group_by(motif) %>%
#   do(RunFisherDHS(., S.collapse, sitecounts.hash))
# 
# # crashes
# # start <- Sys.time()
# # dist.filt <- 1000
# # S.collapse.lf <- subset(S.sub.lf, dist <= dist.filt) %>%
# #   group_by(gene, peak) %>%
# #   do(CollapseDat(., indx = 4, tissue = "Liver", non.tissue = "Flat", flat.style = "normal"))
# # print(Sys.time() - start)
