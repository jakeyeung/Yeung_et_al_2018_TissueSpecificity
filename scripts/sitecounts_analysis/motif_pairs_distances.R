# 2015-01-05
# 

library(parallel)
library(ggplot2)
library(dplyr)
library(hash)
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")


# Function ----------------------------------------------------------------




# Load --------------------------------------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)
load("Robjs/S.collapse.liver.dist1000.Robj", verbose=T)
load("Robjs/N.long.promoters_500.Robj", verbose=T)
# load("Robjs/N.long.all_genes.all_signif_motifs.tmp.Robj", verbose=T)
load("Robjs/N.long.all_genes.all_signif_motifs.Robj", verbose=T)
load("Robjs/S.collapse.liverflat.liver_peaks_only.Robj", verbose=T)  # liver peaks for flat genes
# load("Robjs/N.long.liver_genes.all_motifs.100000.Robj", verbose=T)


# Get sitecounts for all motifs -------------------------------------------
# preprocess data so we only load sitecounts > 0.5 into ram


# Distance between ONECUT1 and other motifs -------------------------------

S.livpeaks <- GetOrderedPeaks(subset(S.collapse, peak.type == "Liver"), fits.best)


# Do Adr flat and rhythmic ------------------------------------------------

models.flat <- ""
models.rhyth <- as.character(subset(fits.best, n.rhyth >= 8)$model)
  
adrmodel <- "Adr"

cutoffs <- seq(from = 0.6, to = 0.6, by = 0.2)
N.adrvflat <- RunFisherOnPromoters(N.long, foreground.models = adrmodel, background.models = models.flat, cutoffs = cutoffs)
N.adrvrhyth <- RunFisherOnPromoters(N.long, foreground.models = adrmodel, background.models = models.rhyth, cutoffs = cutoffs)

N.adr <- data.frame(motif = N.adrvflat$motif, flat = N.adrvflat$odds.ratio, rhyth = N.adrvrhyth$odds.ratio)

ggplot(N.adr, aes(x = flat, y = rhyth, label = motif)) + geom_text() + geom_vline(aes(xintercept = 1)) + geom_hline(aes(yintercept = 1)) + 
  xlab("Odds ratio: adr vs bckgrd (bg=promoters of tissue-wide flat genes)") + ylab("Odds ratio: adr vs bckgrd (bg=promoters of tissue-wide rhythmic genes)") +
  theme_bw(24)

# Do Mef2c flat and rhythmic ----------------------------------------------

models.flat <- ""
models.rhyth <- as.character(subset(fits.best, n.rhyth >= 8)$model)

bfatmodel <- "BFAT"

cutoffs <- seq(from = 0.6, to = 0.6, by = 0.2)
N.bfatvflat <- RunFisherOnPromoters(N.long, foreground.models = bfatmodel, background.models = models.flat, cutoffs = cutoffs)
N.bfatvrhyth <- RunFisherOnPromoters(N.long, foreground.models = bfatmodel, background.models = models.rhyth, cutoffs = cutoffs)

N.bfat <- data.frame(motif = N.bfatvflat$motif, flat = N.bfatvflat$odds.ratio, rhyth = N.bfatvrhyth$odds.ratio)

ggplot(N.bfat, aes(x = flat, y = rhyth, label = motif)) + geom_text() + geom_vline(aes(xintercept = 1)) + geom_hline(aes(yintercept = 1)) + 
  xlab("Odds ratio: BFAT vs bckgrd (bg=promoters of tissue-wide flat genes") + ylab("Odds ratio: BFAT vs bckgrd (bg=promoters of tissue-wide rhythmic genes)") +
  theme_bw(24)


# Do Liver for promoters --------------------------------------------------

models.flat <- ""
models.rhyth <- as.character(subset(fits.best, n.rhyth >= 8)$model)

Livermodel <- "Liver"

cutoffs <- seq(from = 0.6, to = 0.6, by = 0.2)
N.Livervflat <- RunFisherOnPromoters(N.long, foreground.models = Livermodel, background.models = models.flat, cutoffs = cutoffs)
N.Livervrhyth <- RunFisherOnPromoters(N.long, foreground.models = Livermodel, background.models = models.rhyth, cutoffs = cutoffs)

N.Liver <- data.frame(motif = N.Livervflat$motif, flat = N.Livervflat$odds.ratio, rhyth = N.Livervrhyth$odds.ratio)

ggplot(N.Liver, aes(x = flat, y = rhyth, label = motif)) + geom_text() + geom_vline(aes(xintercept = 1)) + geom_hline(aes(yintercept = 1)) + 
  xlab("Odds ratio: Liver vs bckgrd (bg=promoters of tissue-wide flat genes") + ylab("Odds ratio: Liver vs bckgrd (bg=promoters of tissue-wide rhythmic genes)") +
  theme_bw(24)


# Motifs in DHS -----------------------------------------------------------

# what pairs with ONECUT1?

liv.genes <- as.character(subset(fits.best, model == "Liver")$gene)
liv.peaks <- S.livpeaks$peak

flat.genes <- as.character(subset(fits.best, model == "")$gene)
flat.peaks <- subset(S.collapse.livflat, peak.type == "Liver" & model == "Flat")$peak

N.sub <- subset(N.long.filt, gene %in% liv.genes & peak %in% liv.peaks)
N.sub.flat <- subset(N.long.filt, gene %in% flat.genes & peak %in% flat.peaks)

N.sub$motif <- factor(N.sub$motif, levels = unique(N.sub$motif))
N.sub$gene <- factor(N.sub$gene, levels = unique(N.sub$gene))
N.sub$peak <- factor(N.sub$peak, levels = unique(N.sub$peak))

N.sub.flat$motif <- factor(N.sub.flat$motif, levels = unique(N.sub.flat$motif))
N.sub.flat$gene <- factor(N.sub.flat$gene, levels = unique(N.sub.flat$gene))
N.sub.flat$peak <- factor(N.sub.flat$peak, levels = unique(N.sub.flat$peak))


jmotif <- "ONECUT1,2.p2"
jmotif <- "ONECUT1,2.p2"

coop.score <- N.sub %>%
  group_by(peak) %>%
  do(CountDistFromMotif(., jmotif))

coop.score.flat <- N.sub.flat %>%
  group_by(peak) %>%
  do(CountDistFromMotif(., jmotif))

# coop.score.gene <- N.sub %>%
#   group_by(gene) %>%
#   do(CountDistFromMotif(., jmotif))



# ONECUT1 and CUX2 --------------------------------------------------------

ggplot(subset(coop.score.flat, motif != jmotif), aes(x = rel.pos)) + geom_histogram(binwidth = 100) + facet_wrap(~motif) + xlim(c(-1000, 1000))
ggplot(subset(coop.score.flat, motif == "RORA.p2"), aes(x = rel.pos)) + geom_histogram(binwidth = 100) + facet_wrap(~motif) + xlim(c(-1000, 1000))

ggplot(data.frame(subset(coop.score.gene, motif == "CUX2.p2")), aes(x = rel.pos)) + geom_histogram(binwidth = 5) + facet_wrap(~motif) + xlim(c(-100, 100))
ggplot(data.frame(subset(coop.score.gene, motif == "RORA.p2")), aes(x = rel.pos)) + geom_histogram(binwidth = 100) + facet_wrap(~motif) + xlim(c(-1000, 1000))
ggplot(subset(coop.score.gene, motif != jmotif), aes(x = rel.pos)) + geom_histogram(binwidth = 100) + facet_wrap(~motif) + xlim(c(-1000, 1000))

ggplot(subset(coop.score, motif == "CUX2.p2"), aes(x = rel.pos)) + geom_histogram(binwidth = 10) + facet_wrap(~motif)
ggplot(subset(coop.score.flat, motif == "FOXA2.p3"), aes(x = rel.pos)) + geom_histogram(binwidth = 10) + facet_wrap(~motif)
ggplot(subset(coop.score, motif == "ONECUT1,2.p2"), aes(x = rel.pos)) + geom_histogram(binwidth = 10) + facet_wrap(~motif)
ggplot(subset(coop.score, motif == "RXRG_dimer.p3"), aes(x = rel.pos)) + geom_histogram(binwidth = 10) + facet_wrap(~motif)
ggplot(subset(coop.score, motif != jmotif & motif != "CUX2.p2"), aes(x = rel.pos)) + geom_histogram(binwidth = 10) + facet_wrap(~motif)


# Combine liver and flat genes to compare distributions -------------------

coop.score$peak.type <- "Liver"
coop.score.flat$peak.type <- "Flat"

coop.score.all <- rbind(coop.score, coop.score.flat)

# histograms
ggplot(subset(coop.score.all, motif == "CUX2.p2"), aes(x = rel.pos, colour = peak.type)) + geom_histogram(binwidth = 1, fill = "white", position = "identity", alpha = 0.5) + facet_wrap(~motif) + xlim(c(-20, 20))
ggplot(subset(coop.score.all, motif == "HNF1A.p2"), aes(x = rel.pos, colour = peak.type)) + geom_histogram(binwidth = 1, fill = "white", position = "identity", alpha = 0.5) + facet_wrap(~motif) + xlim(c(-20, 20))
ggplot(subset(coop.score.all, motif == "FOXA2.p3"), aes(x = rel.pos, colour = peak.type)) + geom_histogram(binwidth = 10, fill = "white", position = "identity", alpha = 0.5) + facet_wrap(~motif) + xlim(c(-100, 100))
ggplot(subset(coop.score.all, motif == "ZBTB16.p2"), aes(x = rel.pos, colour = peak.type)) + geom_histogram(binwidth = 5, fill = "white", position = "identity", alpha = 0.5) + facet_wrap(~motif) + xlim(c(-100, 100))
ggplot(subset(coop.score.all, motif != jmotif), aes(x = rel.pos, colour = peak.type)) + geom_histogram(binwidth = 1, fill = "white", position = "identity", alpha = 0.5) + facet_wrap(~motif) + xlim(c(-20, 20))
ggplot(subset(coop.score.all, motif != jmotif), aes(x = rel.pos, colour = peak.type)) + geom_histogram(binwidth = 1, fill = "white", position = "identity", alpha = 0.5) + facet_wrap(~motif) + xlim(c(-20, 20))
ggplot(subset(coop.score.all, motif != jmotif), aes(x = rel.pos, colour = peak.type)) + geom_histogram(binwidth = 1, fill = "white", position = "identity", alpha = 0.5) + facet_wrap(~motif) + xlim(c(-100, 100)) + ylim(c(0, 50))
# density plots
ggplot(subset(coop.score.all, motif != jmotif), aes(x = rel.pos, colour = peak.type)) + geom_density() + facet_wrap(~motif)
ggplot(subset(coop.score.all, motif == "HNF1A.p2"), aes(x = rel.pos, colour = peak.type)) + geom_density() + facet_wrap(~motif)
ggplot(subset(coop.score.all, motif == "FOXA2.p3"), aes(x = rel.pos, colour = peak.type)) + geom_density() + facet_wrap(~motif)
ggplot(subset(coop.score.all, motif == "HMGA1,2.p2"), aes(x = rel.pos, colour = peak.type)) + geom_density() + facet_wrap(~motif)
ggplot(subset(coop.score.all, motif == "RORA.p2"), aes(x = rel.pos, colour = peak.type)) + geom_density() + facet_wrap(~motif)
ggplot(subset(coop.score.all, motif == "POU1F1.p2"), aes(x = rel.pos, colour = peak.type)) + geom_density() + facet_wrap(~motif)
ggplot(subset(coop.score.all, motif == "CUX2.p2"), aes(x = rel.pos, colour = peak.type)) + geom_density() + facet_wrap(~motif)

