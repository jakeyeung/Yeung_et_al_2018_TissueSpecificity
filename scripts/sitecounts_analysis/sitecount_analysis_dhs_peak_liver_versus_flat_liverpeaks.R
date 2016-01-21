# 2015-12-10
# Jake Yeung
# same as dhs_peak_gene_body.R but across all motifs.

start <- Sys.time()
library(hash)
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("/home/yeung/projects/tissue-specificity")

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
print("Getting liver genes")

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
S.sub <- S.sub %>%
  group_by(peak, tissue) %>%
  summarise(signal.cut = unique(signal.cut),
            gene = paste(gene, collapse = ";"),
            dist = unique(dist),
            signal.cut = unique(signal.cut))

dist.filt <- 1000
S.collapse <- subset(S.sub, dist <= dist.filt) %>%
  group_by(gene, peak) %>%
  do(CollapseDat(., indx = 4, tissue = "Liver", non.tissue = "Flat", flat.style = "normal"))

# Add flat genes ----------------------------------------------------------
print("Adding flat genes")

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


dist.filt <- 1000
S.collapse.flat <- subset(S.sub.flat, dist <= dist.filt) %>%
  group_by(gene, peak) %>%
  do(CollapseDat(., indx = 4, tissue = "Liver", non.tissue = "Flat", flat.style = "normal"))

save(S.collapse.flat, file = "Robjs/S.collapse.flat.dist1000.Robj")

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
  do(RunFisherDHS(., S.collapse, sitecounts.hash, jmodel.col = "model"))

# plot volcano

# RunFisherDHS(subset(N.long.all, motif == "HIF1A.p2"), S.collapse, sitecounts.hash, jmodel.col = "model", jshow.table = TRUE)

pdf("plots/sitecounts_enrichment_liver_vs_flat_liverpeaks_dhs.pdf")
ggplot(tests.motif.livflat, aes(x = odds.ratio, y = -log10(p.value), label = motif)) + geom_text()
dev.off()

print(Sys.time() - start)
