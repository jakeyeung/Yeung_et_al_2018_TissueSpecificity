# 2015-12-08
# Jake Yeung
# site count analysis on gene body to do motif enrichment

library(ggplot2)
library(mixtools)
library(dplyr)
library(hash)
source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")

# Functions ---------------------------------------------------------------

ReadDHSData <- function(path, tissues, cnames, normalize = TRUE, outlong = TRUE){
  if (missing(path)){
    path <- "/home/yeung/data/tissue_specificity/motevo_dhs/dhs_signal/dhs_signal_windows500.chr.sorted.closest.mat"
  }
  if (missing(tissues)){
    tissues <- c("Cere", "Heart", "Kidney", "Liver", "Lung", "Mus")
  }
  if (missing(cnames)){
    cnames <- c("chromo", "start", "end", tissues, "chromo.gene", "start.gene", "end.gene", "gene", "blank", "strand", "dist")
  }
  S <- read.table("/home/yeung/data/tissue_specificity/motevo_dhs/dhs_signal/dhs_signal_windows500.chr.sorted.closest.mat")
  
  colnames(S) <- cnames
  
  if (normalize){
    for (tiss in tissues){
      S[[tiss]] <- 10^6 * S[[tiss]] / sum(S[[tiss]])
    }
  }
  if (outlong){
    signal.vec <- unlist(S[, colnames(S) %in% tissues])
    S <- data.frame(chromo = S$chromo, start = S$start, end = S$end, 
                         peak = paste(paste(S$chromo, S$start, sep = ":"), S$end, sep = "-"), # chr1:7800234-7800734
                         tissue = rep(tissues, each = nrow(S)), 
                         signal = signal.vec, 
                         gene = S$gene, dist = S$dist)
  }
  return(S)
}

ReadSitecountsMotif <- function(path, cnames, show.time = FALSE){
  start <- Sys.time()
  N <- read.table(path)
  if (missing(cnames)){
    cnames <- c("chromo", "start", "end", "motif_peak", "sitecount", "chromo.gene", "start.gene", "end.gene", "gene", "blank", "strand", "dist")
  }
  colnames(N) <- cnames
  
  # split motif_peak into separate columns
  # RORA.p2;mm10_chr1:3001753-3002253 -> RORA.p2
  N$motif <- sapply(N$motif_peak, function(m) strsplit(as.character(m), ";")[[1]][[1]])
  # RORA.p2;mm10_chr1:3001753-3002253 -> chr1:3001753-3002253
  N$peak <- sapply(N$motif_peak, function(m) strsplit(strsplit(as.character(m), ";")[[1]][[2]], "_")[[1]][[2]])
  
  cnames.remove <- c("motif_peak", "chromo.gene", "start.gene", "end.gene", "blank", "strand")
  for (cname in cnames.remove){
    N[[cname]] <- NULL
  }
  if (show.time){
    print(Sys.time() - start)
  }
  return(N)
}

GetCoord <- function(peak, jget = "start"){
  # chr15:85950195-85950695 -> 85950195 or 85950695 depending on "start" or "end"
  # if chromo, return chr15
  if (jget == "start"){
    jgeti <- 1
  } else if(jget == "end"){
    jgeti <- 2
  } else if(jget == "chromo"){
    return(strsplit(peak, ":")[[1]][[1]])
  } else {
    print(paste("jget must be start or end", jget))
  }
  return(as.numeric(strsplit(strsplit(peak, ":")[[1]][[2]], "-")[[1]][[jgeti]]))
}

GetPeakDistance <- function(peak1, peak2){
  # Get distance between two peaks
  # check they are in same chromosomes
  if (GetCoord(peak1, jget = "chromo") != GetCoord(peak2, jget = "chromo")){
    print("Not on same chromosomes")
    return(NA)
  }
  start1 <- GetCoord(peak1, "start")
  start2 <- GetCoord(peak2, "start")
  end1 <- GetCoord(peak1, "end")
  end2 <- GetCoord(peak2, "end")
  
  dist.min <- min((end1 - start2), (start1 - end2))
  # handle negatives
  return(max(dist.min, 0))
}

# Load --------------------------------------------------------------------


start <- Sys.time()
# N <- read.table("data/sitecounts/motevo_by_peaks_dhs_gene_bodies/merged.closest.bed", nrows = 10)  # 30 GB
S <- read.table("/home/yeung/data/tissue_specificity/motevo_dhs/dhs_signal/dhs_signal_windows500.chr.sorted.closest.mat")
print(Sys.time() - start)

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

FindCutoffLong <- function(dat, signal.col = "signal", jlambdas = c(0.8, 0.2), jmus = c(-1.5, 0.8), 
                           log2.trans = TRUE, pseudo = 1e-2, take.frac = 1, jshow.fig = FALSE){
  # take.frac: sample to a fraction of the data
  if (take.frac < 1){
    dat <- dat[sample(length(dat[[signal.col]]), size = take.frac * nrow(dat), replace = F), ]
  } else if (take.frac > 1){
    print(paste("take.frac must be less than or equal to 1"))
  }
  if (log2.trans){
    cutoff <- FindCutoff(x = log2(dat[[signal.col]] + pseudo), lambdas = jlambdas, mus = jmus, k = 2, show.fig = jshow.fig)
  } else {
    cutoff <- FindCutoff(x = dat[[signal.col]], lambdas = jlambdas, mus = jmus, k = 2, show.fig = jshow.fig)
  }
  if (log2.trans){
    # return in normal scale
    return(data.frame(cutoff = 2^cutoff$maximum))
  } else {
    return(data.frame(cutoff = cutoff$maximum))
  }
}

test <- subset(S.long, tissue == "Heart")
jcut <- FindCutoffLong(test, jlambdas = c(0.7, 0.3), jmus = c(-4, 0), take.frac = 0.001, jshow.fig = TRUE)

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

CollapseDat <- function(dat, tissue, non.tissue = "Flat"){
  # collapse dat into either "tissue" or "non-tissue"
  indx <- which(dat$tissue == tissue)
  tiss.sig <- dat$signal.cut[indx]
  others.sig <- dat$signal.cut[-indx]  # vec
  # check signal.cut has 1 in tissue and 0 in all others
  if (tiss.sig == 1 & max(others.sig) == 0){
    return(data.frame(peak.type = tissue))
  } else if (tiss.sig == 0 & max(others.sig) == 1){
    return(data.frame(peak.type = non.tissue))
  } else {
    return(data.frame())
  }
}

# collapse into liver vs non-liver peaks
S.collapse <- S.sub %>%
  group_by(peak) %>%
  do(CollapseDat(., tissue = "Liver", non.tissue = "Flat"))

# add sitecount info 
N.RORA.sub <- subset(N.RORA, gene %in% liver.genes)
sitecounts.hash <- hash(as.character(N.RORA.sub$peak), N.RORA.sub$sitecount)

S.collapse$sitecount <- sapply(S.collapse$peak, function(p){
  s <- sitecounts.hash[[p]]
  if (is.null(s)){
    return(0)
  } else {
    return(s)
  }
}

FisherTestSitecounts(dat = N.ROR)

