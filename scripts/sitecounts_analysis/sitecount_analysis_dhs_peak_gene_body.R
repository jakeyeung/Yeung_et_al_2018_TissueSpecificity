# 2015-12-08
# Jake Yeung
# site count analysis on gene body to do motif enrichment

library(ggplot2)
library(mixtools)
library(dplyr)
library(hash)
source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
source("scripts/functions/SitecountsFunctions.R")

# Functions ---------------------------------------------------------------

# Load --------------------------------------------------------------------


start <- Sys.time()
# N <- read.table("data/sitecounts/motevo_by_peaks_dhs_gene_bodies/merged.closest.bed", nrows = 10)  # 30 GB
S <- read.table("/home/yeung/data/tissue_specificity/motevo_dhs/dhs_signal/dhs_signal_windows500.chr.sorted.closest.mat")


load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")

sitecounts.dir <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed"
sitecounts.path <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed/RORA.closest.bed"
N.RORA <- ReadSitecountsMotif(sitecounts.path)
N.ONECUT <- ReadSitecountsMotif(file.path(sitecounts.dir, "ONECUT1,2.closest.bed"), show.time = TRUE)
# N.RORA <- read.table(sitecounts.path)

# Add colnames ------------------------------------------------------------

tissues <- c("Cere", "Heart", "Kidney", "Liver", "Lung", "Mus")
cnames <- c("chromo", "start", "end", tissues, "chromo.gene", "start.gene", "end.gene", "gene", "blank", "strand", "dist")
colnames(S) <- cnames


# Normalized dat ----------------------------------------------------------

for (tiss in tissues){
  S[[tiss]] <- 10^6 * S[[tiss]] / sum(S[[tiss]])
}

signal.vec <- unlist(S[, colnames(S) %in% tissues])
S.long <- data.frame(chromo = S$chromo, start = S$start, end = S$end, 
                     peak = paste(paste(S$chromo, S$start, sep = ":"), S$end, sep = "-"), # chr1:7800234-7800734
                     tissue = rep(tissues, each = nrow(S)), 
                     signal = signal.vec, 
                     gene = S$gene, dist = S$dist)


# Regions of interest -----------------------------------------------------

liver.genes <- subset(fits.best, model == "Liver")$gene

S.sub <- subset(S.long, gene %in% liver.genes)

# check for bias across tissues
ggplot(S.long[sample(x = 1:nrow(S.long), size = 0.01 * nrow(S.long)), ], aes(x = log10(signal))) + geom_density() + facet_wrap(~tissue)

# lets look at Celsr1 near chr15:85,959,964-85,964,941
jgene <- "Ube2u"
jchromo <- "chr15"
startmin <- 85959964
endmin <- 85964941

S.subsub <- subset(S.sub, gene == jgene)
N.sub.onecut <- subset(N.ONECUT, gene == jgene)
# N.sub <- subset(N.RORA, gene == jgene & start > startmin)
# N.sub <- subset(N.RORA, chromo == jchromo & start > startmin - 10000 & end < endmin + 10000)
N.sub <- subset(N.RORA, gene == jgene)

# show peaks
ggplot(S.subsub, aes(xmin = start, xmax = end, ymin = -0.5, ymax = 0.5, alpha = signal)) + geom_rect() + 
  geom_rect(aes(xmin = start, xmax = end, ymin = 0.75, ymax = 1.75, alpha = sitecount), data = N.sub.onecut) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 2, ymax = 3, alpha = sitecount), data = N.sub)

ggplot(N.sub, aes(xmin = start, xmax = end, ymin = -0.5, ymax = 0.5)) + geom_rect() + geom_rect(aes(xmin = start, xmax = end, ymin = 0.75, ymax = 1.75), data = N.sub.onecut)
# show peaks 


# Enrichment time ---------------------------------------------------------

# Find cutoff
pseudo <- 1e-2
cutoff <- -2
ggplot(S.long[sample(x = 1:nrow(S.long), size = 0.01 * nrow(S.long)), ], aes(x = log2(signal + pseudo))) + geom_density() + facet_wrap(~tissue) + geom_vline(xintercept = cutoff)

test <- subset(S.long, tissue == "Heart")
jcut <- FindCutoffLong(test, jlambdas = c(0.7, 0.3), jmus = c(-4, 0), take.frac = 0.001, jshow.fig = TRUE)

print(Sys.time() - start)

# needs to be tissue-specific probably, try mixtools
S.tissuecutoff <- S.long %>%
  group_by(tissue) %>%
  do(FindCutoffLong(., jlambdas = c(0.7, 0.3), jmus = c(-4, 0), take.frac = 0.003))

# cool
ggplot(S.long[sample(x = 1:nrow(S.long), size = 0.01 * nrow(S.long)), ], aes(x = log2(signal + pseudo))) + geom_density() + 
  geom_vline(aes(xintercept = log2(cutoff)), data = S.tissuecutoff) + facet_wrap(~tissue)

# now do cutoffs: set DHS signals to 0 or 1
cutoffs.tiss <- hash(as.character(S.tissuecutoff$tissue), as.numeric(S.tissuecutoff$cutoff))
S.sub$signal.cut <- mapply(function(s, tiss){
  cutoff.tiss <- cutoffs.tiss[[tiss]]
  if (s >= cutoff.tiss){
    return(1)
  } else {
    return(0)
  }
}, S.sub$signal, as.character(S.sub$tissue))


# collapse into liver vs non-liver peaks
start <- Sys.time()
S.collapse <- S.sub %>%
  group_by(gene, peak) %>%
  do(CollapseDat(., indx = 4, tissue = "Liver", non.tissue = "Flat", flat.style = "all"))
print(Sys.time() - start)

# # add sitecount info 
# N.RORA.sub <- subset(N.RORA, gene %in% liver.genes)
# N.ONECUT.sub <- subset(N.ONECUT, gene %in% liver.genes)

# collapse readcounts for RORA ONECUT 
N.RORA.sub <- subset(N.RORA, gene %in% liver.genes & dist < 1000) %>%
  group_by(motif, peak) %>%
  summarise(sitecount = sum(sitecount))

N.ONECUT.sub <- subset(N.ONECUT, gene %in% liver.genes & dist  < 1000) %>%
  group_by(motif, peak) %>%
  summarise(sitecount = sum(sitecount))

sitecounts.hash <- hash(as.character(N.RORA.sub$peak), N.RORA.sub$sitecount)
sitecounts.onecut.hash <- hash(as.character(N.ONECUT.sub$peak), N.ONECUT.sub$sitecount)

S.collapse$sitecount.rora <- sapply(S.collapse$peak, AssignSitecount, sitecounts.hash)
S.collapse$sitecount.onecut <- sapply(S.collapse$peak, AssignSitecount, sitecounts.onecut.hash)

FisherTestSitecounts(dat = S.collapse, cutoff = 0.5, sitecount.col = "sitecount.rora", model.col = "peak.type", show.table=TRUE)
FisherTestSitecounts(dat = S.collapse, cutoff = 0.5, sitecount.col = "sitecount.onecut", model.col = "peak.type", show.table=TRUE)


# Do I get RORA if I collapse the peaks assigned to a gene? ---------------

flat.peaks <- subset(S.collapse, peak.type == "Flat")$peak
liver.peaks <- subset(S.collapse, peak.type == "Liver")$peak

N.RORA.flat.gene <- subset(N.RORA, gene %in% liver.genes & dist < 1000 & peak %in% flat.peaks) %>%
  group_by(motif, gene) %>%
  summarise(sitecount = sum(sitecount))
N.RORA.flat.gene$peak.type <- "Flat"

N.RORA.liver.gene <- subset(N.RORA, gene %in% liver.genes & dist < 1000 & peak %in% liver.peaks) %>%
  group_by(motif, gene) %>%
  summarise(sitecount = sum(sitecount))
N.RORA.liver.gene$peak.type <- "Liver"

N.RORA.gene <- rbind(N.RORA.flat.gene, N.RORA.liver.gene)


N.ONECUT.flat.gene <- subset(N.ONECUT, gene %in% liver.genes & dist < 1000 & peak %in% flat.peaks) %>%
  group_by(motif, gene) %>%
  summarise(sitecount = sum(sitecount))
N.ONECUT.flat.gene$peak.type <- "Flat"

N.ONECUT.liver.gene <- subset(N.ONECUT, gene %in% liver.genes & dist < 1000 & peak %in% liver.peaks) %>%
  group_by(motif, gene) %>%
  summarise(sitecount = sum(sitecount))
N.ONECUT.liver.gene$peak.type <- "Liver"

N.ONECUT.gene <- rbind(N.ONECUT.flat.gene, N.ONECUT.liver.gene)

FisherTestSitecounts(dat = N.ONECUT.gene, cutoff = 1, sitecount.col = "sitecount", model.col = "peak.type", show.table=TRUE)
