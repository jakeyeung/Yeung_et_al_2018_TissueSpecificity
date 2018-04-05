# 2016-02-19
# get_tissue_spec_peaks.R
# Find tissue specific peaks within liver-specific rhythmic genes
# but compare KL vs Others (use rhythmic genes always)

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

library(hash)
library(dplyr)
library(reshape2)
library(ggplot2)

# Functions ---------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/DataHandlingFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")

FilterKidneyLiverGenes <- function(param.vec, amp.min){
  # return TRUE if Kidney and Liver amplitudes are above amp.min
  liv.amp <- 0
  kid.amp <- 0  # inits
  livkid.amp <- 0
  
  liv.indx <- which(names(param.vec) == "Liver.amp")
  kid.indx <- which(names(param.vec) == "Kidney.amp")
  livkid.indx <- which(names(param.vec) == "Kidney,Liver.amp")
  if (length(liv.indx) == 1){
    liv.amp <- param.vec[liv.indx]
  }
  if (length(kid.indx) == 1){
    kid.amp <- param.vec[kid.indx]
  } 
  if (length(livkid.indx) == 1){
    livkid.amp <- param.vec[livkid.indx]
  }   
  if ((liv.amp > amp.min & kid.amp > amp.min) | livkid.amp > amp.min) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

FilterLiverGenes <- function(param.vec, amp.min){
  # return TRUE if Kidney and Liver amplitudes are above amp.min
  liv.amp <- 0
  
  liv.indx <- which(names(param.vec) == "Liver.amp")
  if (length(liv.indx) == 1){
    liv.amp <- param.vec[liv.indx]
  }
  if (liv.amp > amp.min) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/S.long.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)
load("Robjs/N.long.all_genes.all_signif_motifs.Robj", v=T)

# Get liver genes ---------------------------------------------------------

# find models that are high in either kidney OR liver but not in "tissue-wide" models
# kidliv.genes <- fits.best$gene[sapply(fits.best$param.list, FilterKidneyLiverGenes, amp.min = 0.25)]
kidliv.genes <- fits.best$gene[sapply(fits.best$param.list, FilterKidneyLiverGenes, amp.min = 0.25)]
# remove tissue-wide genes
kidliv.genes <- subset(fits.best, gene %in% kidliv.genes & n.rhyth < 7)$gene
# plot by random 
print(length(kidliv.genes))
PlotGeneAcrossTissues(subset(dat.long, gene == as.character(sample(kidliv.genes, 1))))

jmodels <- as.character(subset(fits.best, gene %in% kidliv.genes)$model)
model.flat <- ""

compare.vec <- c(F, F, T, T, F, F)  # Cere, Heart, Kidney, Liver, Lung, Mus
# jmodels <-c("Liver")
# compare.vec <- c(F, F, F, T, F, F)  # Cere, Heart, Kidney, Liver, Lung, Mus
jgenes <- unique(as.character(subset(fits.best, model %in% jmodels)$gene))
jgenes.bg <- unique(as.character(subset(fits.best, model == model.flat)$gene))
compare.vec.other <- !compare.vec

# shift by 1 log2 unit the fit looks better
log.shift <- 2.5
# log.shift <- 4.5
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
S.sub.collapse <- subset(S.long, gene %in% c(jgenes, jgenes.bg))

# apply is faster than mapply
S.sub.collapse$is.upper <- apply(S.sub.collapse, 1, function(row) IsSignalUpper(as.character(row[5]), as.numeric(row[6]), cutoffs.tiss.upper))
S.sub.collapse$is.lower <- apply(S.sub.collapse, 1, function(row) IsSignalLower(as.character(row[5]), as.numeric(row[6]), cutoffs.tiss.lower))

# all peaks with high signal
S.sub.collapse.peaks <- subset(S.sub.collapse, is.upper == TRUE)
S.sub.collapse.peaks.lower <- subset(S.sub.collapse, is.lower == TRUE)


# get livkid peaks within foreground genes and other peaks from foreground genes
livkid.peaks <- unique(as.character(subset(S.sub.collapse.peaks, gene %in% jgenes & tissue %in% c("Kidney", "Liver") & is.upper == TRUE)$peak))
bg.peaks <- unique(as.character(subset(S.sub.collapse.peaks, gene %in% jgenes & ! tissue %in% c("Kidney", "Liver") & is.upper == TRUE)$peak))
flat.peaks <- unique(as.character(subset(S.sub.collapse.peaks, gene %in% jgenes.bg & tissue %in% c("Kidney", "Liver") & is.upper == TRUE)$peak))

S.sub.collapse.livkid <- subset(S.sub.collapse, peak %in% livkid.peaks) %>%
  group_by(peak, gene) %>%
  summarise(is.tissspec = IsTissSpec(is.upper, compare.vec)) %>%
  filter(is.tissspec == TRUE)
#   summarise(is.tissspec = IsTissSpec(is.upper, compare.vec),
#             is.not.tissspec = IsTissSpec(is.lower, !compare.vec)) %>%
#   filter(is.tissspec == TRUE & is.not.tissspec)

# compare.vec.other <- c(T, T, T, T, T, T)
# compare.vec.other <- c(F, F, F, F, F, F)
compare.vec.other <- c(F, F, F, F, F, F)
S.sub.collapse.bg <- subset(S.sub.collapse, peak %in% bg.peaks) %>%
  group_by(peak, gene) %>%
  summarise(is.tissspec = IsTissSpec(is.upper, compare.vec.other, reverse = TRUE)) %>%
  filter(is.tissspec == TRUE)

S.sub.collapse.flat <- subset(S.sub.collapse, peak %in% flat.peaks) %>%
  group_by(peak, gene) %>%
  summarise(is.tissspec = IsTissSpec(is.upper, compare.vec)) %>%
  filter(is.tissspec == TRUE)

print(paste("N foreground peaks:", length(unique(as.character(S.sub.collapse.livkid$peak)))))
print(paste("N background peaks:", length(unique(as.character(S.sub.collapse.bg$peak)))))
print(paste("N flat peaks:", length(unique(as.character(S.sub.collapse.flat$peak)))))



# Sanity check my peeks ---------------------------------------------------

# signal in liver
jsub.liv <- subset(S.sub.collapse, peak %in% as.character(S.sub.collapse.livkid$peak) & tissue == "Liver")
hist(jsub.liv$signal, breaks = 100)

jsub.kid <- subset(S.sub.collapse, peak %in% as.character(S.sub.collapse.livkid$peak) & tissue == "Kidney")
hist(jsub.kid$signal, breaks = 100)

# # cutoff here
# sig.cut <- 1.881585
# fg.peaks.filt <- as.character(subset(jsub.liv, signal > sig.cut)$peak)
# print(length(fg.peaks.filt))
# S.sub.collapse.livkid <- subset(S.sub.collapse.livkid, peak %in% fg.peaks.filt)
# 

# down sample to about the same size
set.seed(0)
bg.tissspec.peaks <- unique(as.character(S.sub.collapse.bg$peak))
fg.tissspec.peaks <- unique(as.character(S.sub.collapse.livkid$peak))
bg.tissspec.peaks.samp <- sample(x = unique(as.character(S.sub.collapse.bg$peak)), size = length(fg.tissspec.peaks))
S.sub.collapse.bg <- subset(S.sub.collapse.bg, peak %in% bg.tissspec.peaks.samp)

# Make sitecounts ---------------------------------------------------------

head(N.long.filt)
N.sub.livkid <- subset(N.long.filt, peak %in% S.sub.collapse.livkid$peak)
N.sub.other <- subset(N.long.filt, peak %in% S.sub.collapse.bg$peak)
N.sub.flat <- subset(N.long.filt, peak %in% S.sub.collapse.flat$peak)

# run LDA
mat.fg <- dcast(subset(N.sub.livkid), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bg <- dcast(subset(N.sub.other), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.flat <- dcast(subset(N.sub.flat), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)

# bind rows create labels
mat.fgbg <- bind_rows(mat.fg, mat.bg)
mat.fgbg[is.na(mat.fgbg)] <- 0
rownames(mat.fgbg) <- paste(mat.fgbg$peak, mat.fgbg$gene, sep = ";"); mat.fgbg$peak <- NULL; mat.fgbg$gene <- NULL
labels <- c(rep(1, nrow(mat.fg)), rep(2, nrow(mat.bg)))

out <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  

BoxplotLdaOut(out, jtitle = "Single factor separation")
PlotLdaOut(out)
m.singles <- SortLda(out)

# compare with flat
mat.fgflat <- bind_rows(mat.fg, mat.flat)
mat.fgflat[is.na(mat.fgflat)] <- 0
rownames(mat.fgflat) <- make.unique(paste(mat.fgflat$peak, mat.fgflat$gene, sep = ";")); mat.fgflat$peak <- NULL; mat.fgflat$gene <- NULL
labels.flat <- c(rep(1, nrow(mat.fg)), rep(2, nrow(mat.flat)))
out.flat <- PenalizedLDA(mat.fgflat, labels.flat, lambda = 0.1, K = 1, standardized = FALSE)  
BoxplotLdaOut(out.flat, jtitle = "Single factor separation: vs flat")
PlotLdaOut(out.flat)
m.fgflat <- SortLda(out.flat)

# do crosses
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)
mat.fgbg.cross <- CrossProduct(mat.fgbg, remove.duplicates = TRUE)
dim(mat.fgbg.cross)
# add single factors
mat.fgbg.cross <- cbind(mat.fgbg, mat.fgbg.cross)
# remove columns with 0 variance 
mat.fgbg.cross[which(colSums(mat.fgbg.cross) == 0)] <- list(NULL)

# jlambda <- 0.029  # kidliv
jlambda <- 0.015  # liv only
out.cross <- PenalizedLDA(mat.fgbg.cross, labels, lambda = jlambda, K = 1, standardized = FALSE)
m <- SortLda(out.cross)
print(length(m))
BoxplotLdaOut(out.cross, jtitle = "Cross product separation")
PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE)

# do crosses with flat
colnames(mat.fgflat) <- sapply(colnames(mat.fgflat), RemoveP2Name)
mat.fgflat.cross <- CrossProduct(mat.fgflat, remove.duplicates=TRUE)
dim(mat.fgflat.cross)
mat.fgflat.cross <- cbind(mat.fgflat, mat.fgflat.cross)
# remove columns with 0 variance 
mat.fgflat.cross[which(colSums(mat.fgflat.cross) == 0)] <- list(NULL)

jlambda <- 0.015  # liv only
out.flat.cross <- PenalizedLDA(mat.fgflat.cross, labels.flat, lambda = jlambda, K = 1, standardized = FALSE)
m.flat.cross <- SortLda(out.flat.cross)
print(length(m.flat.cross))
BoxplotLdaOut(out.flat.cross, jtitle = "Cross product separation")
PlotLdaOut(out.flat.cross, take.n = 50, from.bottom = TRUE)

# do across 3 labels
mat.fgbgflat <- bind_rows(mat.fg, mat.bg, mat.flat)
mat.fgbgflat[is.na(mat.fgbgflat)] <- 0
rownames(mat.fgbgflat) <- make.unique(paste(mat.fgbgflat$peak, mat.fgbgflat$gene, sep = ";")); mat.fgbgflat$peak <- NULL; mat.fgbgflat$gene <- NULL
labels.bgflat <- c(rep(1, nrow(mat.fg)), rep(2, nrow(mat.bg)), rep(3, nrow(mat.flat)))
out.bgflat <- PenalizedLDA(mat.fgbgflat, labels.bgflat, lambda = 0.1, K = 2, standardized = FALSE)  
# BoxplotLdaOut(out.bgflat, jtitle = "Single factor separation: vs flat and same")
PlotLdaOut2D(out.bgflat, jdim = 2)
m.fgbgflat <- SortLda(out.bgflat)

# do 2 labels but bg and flat are merged
mat.fgbgflat.merged <- bind_rows(mat.fg, mat.bg, mat.flat)
mat.fgbgflat.merged[is.na(mat.fgbgflat.merged)] <- 0
rownames(mat.fgbgflat.merged) <- make.unique(paste(mat.fgbgflat.merged$peak, mat.fgbgflat.merged$gene, sep = ";")); mat.fgbgflat.merged$peak <- NULL; mat.fgbgflat.merged$gene <- NULL
labels.bgflat.merged <- c(rep(1, nrow(mat.fg)), rep(2, (nrow(mat.bg) + nrow(mat.flat))))
out.bgflat.merged <- PenalizedLDA(mat.fgbgflat.merged, labels.bgflat.merged, lambda = 0.1, K = 1, standardized = FALSE)  
BoxplotLdaOut(out.bgflat.merged, jtitle = "Single factor separation: vs flat and same")
PlotLdaOut(out.bgflat.merged)
m.fgbgflat.merged <- SortLda(out.bgflat.merged)


# Are you sure this is real? ----------------------------------------------

# how many DBP?
mat.fg.sums <- apply(mat.fg[, 3:ncol(mat.fg)], 2, sum)
mat.bg.sums <- apply(mat.bg[, 3:ncol(mat.bg)], 2, sum)
length(mat.fg.sums)

plot(mat.fg.sums, pch=".")
text(mat.fg.sums, labels = names(mat.fg.sums))

plot(mat.bg.sums, pch=".")
text(mat.bg.sums, labels = names(mat.bg.sums))

# and for crosses (take subset)
mat.fg.cross <- mat.fgbg.cross[1:nrow(mat.fg), ]
mat.bg.cross <- mat.fgbg.cross[(nrow(mat.fg)+1):(nrow(mat.fg)+nrow(mat.bg)), ]
mat.flat.cross <- mat.fgflat.cross[(nrow(mat.fg)+1):(nrow(mat.fg)+nrow(mat.flat)), ]
mat.fg.cross.sums <- apply(mat.fg.cross[, 3:ncol(mat.fg.cross)], 2, sum)
mat.bg.cross.sums <- apply(mat.bg.cross[, 3:ncol(mat.bg.cross)], 2, sum)
mat.flat.cross.sums <- apply(mat.flat.cross[, 3:ncol(mat.flat.cross)], 2, sum)


# take top 200
mat.fg.cross.sums.filt <- head(sort(mat.fg.cross.sums, decreasing = TRUE), n = 200)

# how many DBP-HNF4A hits?
jmotif <- "DBP;HNF4A_NR2F1,2"
jmotif <- "RORA;HNF4A_NR2F1,2"
jmotif <- "HNF4A_NR2F1,2;RORA"
jmotif <- "HNF4A_NR2F1,2;RXRG_dimer"
jmotif <- "HNF4A_NR2F1,2"
jmotif <- "ADNP_IRX_SIX_ZHX;DBP"
jmotif <- "HNF4A_NR2F1,2"
jmotif <- "DBP"
jmotif <- "RORA"
jmotif <- "IKZF1;RORA"
fg.hits <- rownames(mat.fg.cross)[which(mat.fg.cross[, c(jmotif)] > 0)]
bg.hits <- rownames(mat.bg.cross)[which(mat.bg.cross[, c(jmotif)] > 0)]
flat.hits <- rownames(mat.flat.cross)[which(mat.flat.cross[, c(jmotif)] > 0)]
print(paste("hits of", jmotif, "FG", length(fg.hits)))
print(paste("hits of", jmotif, "BG", length(bg.hits)))
print(paste("hits of", jmotif, "Flat", length(flat.hits)))

# which pair has most hits?
n.hits <- apply(mat.fg.cross[, grepl(";", colnames(mat.fg.cross))], 2, function(jcol) length(which(jcol > 0)))
print(head(sort(n.hits, decreasing = TRUE), n = 20))

# top hits with RORA
head(sort(n.hits[grepl("RORA", names(n.hits))], decreasing=TRUE))

# Write peaks to output  --------------------------------------------------

# source("scripts/functions/PlotUCSC.R")  # CoordToBed
# 
# fg.peaks.final <- unique(as.character(S.sub.collapse.livkid$peak))
# fg.bed <- lapply(fg.peaks.final, CoordToBed)
# fg.bed <- do.call(rbind, fg.bed)
# fg.bed$blank <- "."
# fg.bed$id <- paste0("mm10_", fg.peaks.final)
# write.table(fg.bed, file = "bedfiles/lda_analysis/fg_kidliv_more_genes_peaks.tmp.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# 
# bg.peaks.final <- unique(as.character(S.sub.collapse.bg$peak))
# bg.bed <- lapply(bg.peaks.final, CoordToBed)
# bg.bed <- do.call(rbind, bg.bed)
# bg.bed$blank <- "."
# bg.bed$id <- paste0("mm10_", bg.peaks.final)
# write.table(fg.bed, file = "bedfiles/lda_analysis/bg_kidliv_more_genes_peaks.tmp.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# 
# # e.g.,
# jcoord <- "chr7:99971744-99972244"  # Rnf169
# PlotGeneAcrossTissues(subset(dat.long, gene=="Rnf169"))

