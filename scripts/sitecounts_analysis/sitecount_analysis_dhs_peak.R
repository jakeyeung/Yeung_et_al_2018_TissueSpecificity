# 2015-11-30
# Jake Yeung
# sitecount_analysis_dhs_peaks.R

library(dplyr)
library(ggplot2)
library(hash)

source("scripts/functions/FisherTestSitecounts.R")

AddAnnots <- function(mat){
  # chr1:3204628-3205128;Xkr4 extract chromo, start, end, gene
  for (func in c(AddChromo, AddStart, AddEnd, AddGene)){
    mat <- func(mat)
  }
  # add chr1:3204628-3205128
  mat$id <- sapply(rownames(mat), function(r) strsplit(r, ";")[[1]][[1]])
  return(mat)
}

AddChromo <- function(mat){
  mat$chromo <- sapply(rownames(mat), function(r) strsplit(r, ":")[[1]][[1]])
  return(mat)
}
AddStart <- function(mat){
  mat$start <- sapply(rownames(mat), function(r) strsplit(strsplit(r, ":")[[1]][[2]], "-")[[1]][[1]])
  return(mat)
}
AddEnd <- function(mat){
  mat$end <- sapply(rownames(mat), function(r) strsplit(strsplit(r, "-")[[1]][[2]], ";")[[1]][[1]])
  return(mat)
}
AddGene <- function(mat){
  mat$gene <- sapply(rownames(mat), function(r) strsplit(r, ";")[[1]][[2]])
  return(mat)
}

MakeSubFromGenelists <- function(dat, genelists, genelist_names){
  # Make sub from genelists and add column "model" describing name of genelist
  if (length(genelists) != length(genelist_names)) warning("Length of genelists and length of names must be equal")
  
  dat.new <- data.frame()  # init
  for (i in 1:length(genelists)){
    dat.sub <- subset(dat, gene %in% genelists[[i]])
    dat.sub$model <- genelist_names[[i]]
    dat.new <- rbind(dat.new, dat.sub)
  }
  return(dat.new)
}

IsSpecific <- function(dat, cutoff){
  # Check if dat chunk (peak between Liver and Other) are 
  # specific to Liver or Other, if so then return with
  # the tissue that it is rhythmic
  if (range(dat$dhs_signal)[1] < cutoff & range(dat$dhs_signal)[2] > cutoff){
    jtiss <- dat$tissue[which.max(dat$dhs_signal)]
  } else {
    jtiss <- NA
  }
  if (is.na(jtiss)){
    return(data.frame())
  } else {
    return(data.frame(tissue = jtiss, id = dat$id[1], gene = dat$gene[1], dhs_signal.max = max(dat$dhs_signal)))
  }
}

# Load --------------------------------------------------------------------

N <- read.table("/home/yeung/projects/tissue-specificity/data/sitecounts/motevo_by_peaks_dhs/sitecounts.filtered.cutoff0.255.mat")
dhs <- read.table("/home/yeung/projects/tissue-specificity/data/sitecounts/motevo_by_peaks_dhs/dhs_signal.filtered.cutoff0.255.mat")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)

N.backup <- N

N <- AddAnnots(N)
dhs <- AddAnnots(dhs)

genes.liv <- as.character(subset(fits.best, model == "Liver")$gene)
genes.flat <- as.character(subset(fits.best, model == "")$gene)

N.sub <- MakeSubFromGenelists(N, list(genes.liv, genes.flat), list("Liver", "Flat"))
dhs.sub <- MakeSubFromGenelists(dhs, list(genes.liv, genes.flat), list("Liver", "Flat"))

# take subset will make life easier i think
N.long <- melt(N.sub, value.name = "sitecount", id.vars = c("chromo", "start", "end", "gene", "model", "id"), variable.name = "motif")
dhs.long <- melt(dhs.sub, value.name = "dhs_signal", id.vars = c("chromo", "start", "end", "gene", "model", "id"), variable.name = "tissue")

dhs.long <- subset(dhs.long, tissue %in% c("Cerebellum", "Heart", "Kidney", "Liver", "Lung", "SkeletalMuscle"))

# Compress dhs.long to Liver and Other
dhs.long$tissue2 <- as.factor(sapply(dhs.long$tissue, function(c){
  if (c == "Liver"){
    return("Liver")
  } else {
    return("Other")
  }
}))

dhs.long$id <- as.factor(dhs.long$id)

dhs.long2 <- dhs.long %>%
  group_by(tissue2, id, gene) %>%
  summarise(dhs_signal = max(dhs_signal))
dhs.long2$tissue <- dhs.long2$tissue2; dhs.long2$tissue2 <- NULL

# use cutoff log10(dhs_signal) = -0.9
cutoff <- 10^-0.9

# filter for dhs sites that are specific to either Liver or Other
dhs.long2 <- dhs.long2 %>%
  group_by(id, gene) %>%
  do(IsSpecific(., cutoff))

# get only tissue-specific peaks
ids.tissuespec <- dhs.long2$id

N.long.sub <- subset(N.long, id %in% ids.tissuespec)

# label each peak by whether it is specific to liver or other
tissue.hash <- hash(as.character(dhs.long2$id), as.character(dhs.long2$tissue))

N.long.sub$peak.tiss <- as.factor(sapply(as.character(N.long.sub$id), function(id) tissue.hash[[id]]))

jmotif <- "ONECUT1.2.p2"
jmotif <- "RORA.p2"
FisherTestSitecounts(dat = subset(N.long.sub, motif == jmotif), cutoff = 0.5, sitecount.col = "sitecount", model.col = "peak.tiss", show.table = TRUE)

jgene <- "Celsr1"  # peak at chr15:85,961,713-85,962,101
jmotif <- "RORA.p2"
test <- subset(N.long.sub, gene == "Celsr1" & motif == "RORA.p2")
test.full <- subset(N.long, gene == "Celsr1" & motif == "RORA.p2")

test.sum <- test %>%
  group_by(peak.tiss) %>%
  summarise(sitecount.mean = mean(sitecount))



ggplot(test, aes(x = peak.tiss, y = sitecount)) + geom_boxplot()

# merge
# Ndhs.long <- merge(N.long, dhs.long)  # too long dont do

# Trying to do analyses on merged DHS profiles are no good!
# 
# # N <- read.table("/home/yeung/projects/tissue-specificity/data/Asitecounts/motevo_by_peaks_by_tissue/all_sites.by_tissue.closest.bed")  # takes forever!
# # N <- read.table("/home/yeung/projects/tissue-specificity/data/sitecounts/motevo_by_peaks_by_tissue/all_sites.closest.smaller.bed")
# load("Robjs/sitecounts_by_motif2.Robj", verbose=T)
# load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
# 
# # N$motif <- sapply(as.character(N$motif), function(m) strsplit(m, ";")[[1]][[1]])
# # save(N, file = "Robjs/sitecounts_by_motif2.Robj")
# # Find liver-specific peaks -----------------------------------------------
# 
# genes.liv <- as.character(subset(fits.best, model == "Liver")$gene)
# genes.flat <- as.character(subset(fits.best, model == "")$gene)
# 
# N.sub <- subset(N, gene %in% c(genes.liv, genes.flat))
# 
# 
# # Motif enrichment --------------------------------------------------------
# 
# jtiss <- "Liver"
# 
# N.sub.liv <- subset(N.sub, gene %in% genes.liv)
# N.sub.flat <- subset(N.sub, gene %in% genes.flat)
# 
# N.sub.liv.liv <- subset(N.sub.liv, tissue == jtiss)
# N.sub.flat.liv <- subset(N.sub.liv, tissue != jtiss)
# 
# N.sub.liv.liv$model <- "LiverRhyth"
# N.sub.flat.liv$model <- "Flat"
# 
# N.sub.test <- rbind(N.sub.liv.liv, N.sub.flat.liv)
# 
# FisherTestSitecounts(dat = subset(N.sub.test, motif == "HNF1A.p2"), cutoff = 0.7, sitecount.col = "sitecount", show.table = TRUE)
# 
# start <- Sys.time()
# cutoffs <- seq(from = 0.4, to = 0.8, by = 0.1)
# N.ftest.ltw.all <- data.frame()
# for (cutoff in cutoffs){
#   print(cutoff)
#   N.ltw.ftest <- N.sub.test %>%
#     group_by(motif) %>%
#     do(FisherTestSitecounts(., cutoff))
#   N.ltw.ftest$cutoff <- cutoff
#   N.ftest.ltw.all <- rbind(N.ftest.ltw.all, N.ltw.ftest)
# }
# print(Sys.time() - start)
# N.ftest.ltw.sum <- N.ftest.ltw.all %>%
#   group_by(motif) %>%
#   summarise(odds.ratio = mean(odds.ratio), p.value = mean(p.value))
# ggplot(N.ftest.ltw.sum, aes(y = -log10(p.value), x = odds.ratio, label = motif)) + geom_point() + geom_text()