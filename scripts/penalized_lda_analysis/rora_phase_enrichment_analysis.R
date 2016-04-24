# 2016-04-22
# Subset plda, also ask how many liver-specific peaks there are in different sets?

rm(list=ls())

start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(parallel)
library(hash)

# Functions ---------------------------------------------------------------


# Constnts ----------------------------------------------------------------

distfilt <- 5000
jcutoff <- 1.5

# Load --------------------------------------------------------------------

load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/N.long.livertwflat.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
load("Robjs/wtko.dat.wtko.hog.commongenes.Robj", v=T)
load("Robjs/wtko.fits.all.long.commongenes.Robj", v=T)
load("Robjs/dat.fit.Robj", v=T)

# Get liv genes -----------------------------------------------------------

# Get genes and peaks -----------------------------------------------------

jmodels <- c("Liver")

amp.min <- 0
jgenes <- as.character(subset(fits.best, model %in% jmodels & amp.avg > amp.min)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)

print(paste("Rhythmic genes:", length(jgenes)))
print(paste("Flat genes:", length(jgenes.flat)))


print(paste("Rhythmic genes:", length(jgenes)))
print(paste("Flat genes:", length(jgenes.flat)))

# get peaks near gene
S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
jpeaks.flat <- as.character(unique(S.sub.flat$peak))

print(paste("number of peaks in liver-specific rhythmic genes", length(jpeaks)))
print(paste("number of peaks in flat genes", length(jpeaks.flat)))

N.sub <- subset(N.long.filt, gene %in% jgenes & peak %in% jpeaks)  # should not have any missing peaks
N.sub.flat <- subset(N.long.filt, gene %in% jgenes.flat & peak %in% jpeaks.flat)


# Tissue-specific peaks ---------------------------------------------------

rhyth.tiss <- c("Liver")
# take peaks with Liver signal greater than cutoff
jtiss <- levels(S.sub$tissue)
tiss.i <- which(jtiss %in% rhyth.tiss)
others.i <- which(!jtiss %in% rhyth.tiss)

S.sub.liverpeaks <- S.sub %>%
  group_by(peak) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
  # filter(mean(zscore[tiss.i]) > jcutoff & mean(zscore[others.i]) < jcutoff)
  filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff)

S.sub.nonliverpeaks <- S.sub %>%
  group_by(peak) %>%
  # filter(mean(zscore[tiss.i]) < jcutoff & mean(zscore[others.i]) > jcutoff)  # any tissue needs to be above cutoff
  # filter(mean(zscore[others.i]) > jcutoff)  # avg of all tissues greater than cutoff
  filter(max(zscore[tiss.i] < jcutoff) & max(zscore[others.i]) > jcutoff)  # any bg tissue above cutoff

jtiss.flat <- levels(S.sub.flat$tissue)
if (identical(jtiss, jtiss.flat) == FALSE){
  print("This shouldnt be necessary")
  tiss.i <- which(jtiss %in% rhyth.tiss)
  others.i <- which(jtiss %in% rhyth.tiss)
}

S.sub.flat.liverpeaks <- S.sub.flat %>%
  group_by(peak) %>%
  # filter(mean(zscore[tiss.i]) > jcutoff & mean(zscore[others.i]) < jcutoff)
  filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff)


# Define peaks ------------------------------------------------------------

liver.peaks.fg <- as.character(unique(S.sub.liverpeaks$peak))
nonliver.peaks.bg <- as.character(unique(S.sub.nonliverpeaks$peak))
liver.peaks.bg <- as.character(unique(S.sub.flat.liverpeaks$peak))
intersect.peaks <- intersect(liver.peaks.fg, liver.peaks.bg)
intersect.peaks.nonliver <- intersect(liver.peaks.fg, nonliver.peaks.bg)
if (length(intersect.peaks.nonliver) == 0){
  warning("Warning, there are intersecting peaks between Liver FG and Non Liver BG!")
}

liver.peaks.fg <- liver.peaks.fg[which(!liver.peaks.fg %in% intersect.peaks)]
liver.peaks.bg <- liver.peaks.bg[which(!liver.peaks.bg %in% intersect.peaks)]
#   liver.peaks.
print(paste("N.peaks fg", length(liver.peaks.fg)))
print(paste("N.peaks bg", length(liver.peaks.bg)))
print(paste("N.peaks flat tiss bg", length(nonliver.peaks.bg)))


# How many genes with ROR? ------------------------------------------------

# not peak, but GENES

liverrhyth.rora <- subset(N.sub, motif == "RORA.p2" & peak %in% liver.peaks.fg) %>%
  group_by(gene) %>%
  summarise(max.rora = max(sitecount)) %>%
  mutate(peak.type = "LiverRhyth")

liverflat.rora <- subset(N.sub.flat, motif == "RORA.p2" & peak %in% liver.peaks.bg) %>%
  group_by(gene) %>%
  summarise(max.rora = max(sitecount)) %>%
  mutate(peak.type = "LiverFlat")

nonliver.rora <- subset(N.sub, motif == "RORA.p2" & peak %in% nonliver.peaks.bg) %>%
  group_by(gene) %>%
  summarise(max.rora = max(sitecount)) %>%
  mutate(peak.type = "NonLiver")

rora.all <- do.call(rbind, list(liverrhyth.rora, liverflat.rora, nonliver.rora))

ggplot(rora.all, aes(x = max.rora)) + geom_histogram() + facet_wrap(~peak.type)

ggplot(rora.all, aes(x = max.rora)) + geom_density() + facet_wrap(~peak.type)


# Phase distributions of liver-specific rhythmic genes --------------------

clock.controlled.genes <- as.character(subset(fits.all.long.wtkohog, model %in% c("WT,Liver", "WT;KO;Liver", "WT;Liver", "Liver"))$gene)
jgenes.filt <- intersect(jgenes, clock.controlled.genes)

rora.genes <- as.character(subset(liverrhyth.rora, max.rora > 0.5)$gene)
ggplot(subset(dat.fit, tissue == "Liver" & gene %in% jgenes), aes(x = phase)) + geom_histogram(bins = 40)
ggplot(subset(dat.fit, tissue == "Liver" & gene %in% as.character(liverrhyth.rora$gene)), aes(x = phase)) + geom_histogram(bins = 40)
ggplot(subset(dat.fit, tissue == "Liver" & gene %in% rora.genes), aes(x = phase)) + geom_histogram(bins = 40)
ggplot(subset(dat.fit, tissue == "Liver" & gene %in% jgenes.filt), aes(x = phase)) + geom_histogram(bins = 40)


