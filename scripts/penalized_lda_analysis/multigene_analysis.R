# 2016-04-04
# Jake Yeung
# multigene_analysis.R
# Assign peaks to many genes: how does this change our analysis?

rm(list=ls())
library(dplyr)
library(reshape2)
library(ggplot2)
library(hash)
library(penalizedLDA)


# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/ShannonEntropy.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")


# Function ----------------------------------------------------------------


SetUpMatForLda <- function(mat.fg, mat.bg){
  mat.fgbg <- bind_rows(mat.fg, mat.bg)
  mat.fgbg[is.na(mat.fgbg)] <- 0
  rownames(mat.fgbg) <- paste(mat.fgbg$peak, mat.fgbg$gene, sep = ";"); mat.fgbg$peak <- NULL; mat.fgbg$gene <- NULL
  labels <- c(rep(1, nrow(mat.fg)), rep(2, nrow(mat.bg)))
  return(list(mat.fgbg = mat.fgbg, labels = labels))
}


# Load --------------------------------------------------------------------


load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/N.long.multigene.distfilt.50000.Robj", v=T)
# subset(N.long.filt, 
load("Robjs/S.long.multigene.filt.50000.Robj", v=T)

# Main --------------------------------------------------------------------

# colnames(N.multigene) <- c("chromo", "start", "end", "motif", "sitecount", "chromo.peak", "start.peak", "end.peak", "gene", "dist")
distfilt <- 2500
amp.min <- 0
rhyth.tiss <- c("Liver")

jgenes <- as.character(subset(fits.best, model == "Liver" & amp.avg > amp.min)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)

# get peaks near gene
S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
jpeaks.flat <- as.character(unique(S.sub.flat$peak))

jcutoff <- 1.5  # arbitrary
# ggplot(S.sub, aes(x = zscore)) + geom_histogram(bins = 100) + facet_wrap(~tissue) + theme_bw(24) + geom_vline(aes(xintercept = jcutoff))
# ggplot(S.sub.flat, aes(x = zscore)) + geom_histogram(bins = 100) + facet_wrap(~tissue) + theme_bw(24) + geom_vline(aes(xintercept = jcutoff))


# Subset sitecounts (distfilt from S) -------------------------------------

N.sub <- subset(N.long.filt, peak %in% jpeaks)  # should not have any missing peaks
N.sub.flat <- subset(N.long.filt, peak %in% jpeaks.flat)


# Take only peaks that are "tissue-specific?" -----------------------------

# take peaks with Liver signal greater than cutoff
jtiss <- levels(S.sub$tissue)
tiss.i <- which(jtiss %in% rhyth.tiss)
others.i <- which(!jtiss %in% rhyth.tiss)

S.sub.liverpeaks <- S.sub %>%
  group_by(peak) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
  filter(zscore[tiss.i] > jcutoff & max(zscore[others.i]) < jcutoff)

S.sub.nonliverpeaks <- S.sub %>%
  group_by(peak) %>%
  filter(zscore[tiss.i] < jcutoff & max(zscore[others.i]) > jcutoff)

jtiss.flat <- levels(S.sub.flat$tissue)
if (identical(jtiss, jtiss.flat) == FALSE){
  # shouldnt be necessary
  tiss.i <- which(jtiss == "Liver")
  others.i <- which(jtiss != "Liver")
}

S.sub.flat.liverpeaks <- S.sub.flat %>%
  group_by(peak) %>%
  filter(zscore[tiss.i] > jcutoff & max(zscore[others.i]) < jcutoff)
  

# Do penalized LDA --------------------------------------------------------

liver.peaks.fg <- as.character(unique(S.sub.liverpeaks$peak))
nonliver.peaks.bg <- as.character(unique(S.sub.nonliverpeaks$peak))
liver.peaks.bg <- as.character(unique(S.sub.flat.liverpeaks$peak))
intersect.peaks <- intersect(liver.peaks.fg, liver.peaks.bg)
print(paste("N.peaks rhyth tiss fg", length(liver.peaks.fg)))
print(paste("N.peaks flat tiss bg", length(nonliver.peaks.bg)))
print(paste("N.peaks rhyth tiss bg", length(liver.peaks.bg)))

liver.peaks.fg <- liver.peaks.fg[which(!liver.peaks.fg %in% intersect.peaks)]
liver.peaks.bg <- liver.peaks.bg[which(!liver.peaks.bg %in% intersect.peaks)]
print(paste("N.peaks fg", length(liver.peaks.fg)))
print(paste("N.peaks bg", length(liver.peaks.bg)))


mat.fg <- dcast(subset(N.sub, peak %in% liver.peaks.fg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bgnonliver <- dcast(subset(N.sub, peak %in% nonliver.peaks.bg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bg <- dcast(subset(N.sub.flat, peak %in% liver.peaks.bg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)

# bind rows create labels

mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bgnonliver)

mat.fgbg <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels

out <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  

BoxplotLdaOut(out, jtitle = "Single factor separation")
PlotLdaOut(out, jtitle = "Single factor loadings")

m.singles <- SortLda(out)


# Do cross products -------------------------------------------------------

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
PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE, jtitle = "Cross product loadings")


# How many RORA-partners happen? ------------------------------------------

n.hits <- apply(mat.fgbg.cross[, grepl(";", colnames(mat.fgbg.cross))], 2, function(jcol) length(which(jcol > 0)))
# print(head(sort(n.hits, decreasing = TRUE), n = 20))

# top hits with RORA
rora.pairs <- sort(n.hits[grepl("RORA", names(n.hits))], decreasing=TRUE)
rora.pairs.filt <- rora.pairs[1:30]
textplot(seq(length(rora.pairs.filt)), rora.pairs.filt, words = names(rora.pairs.filt), xlab = "Index", ylab = "Number of pairs in foreground", cex = 0.9)

onecut.pairs <- sort(n.hits[grepl("ONECUT", names(n.hits))], decreasing=TRUE)
onecut.pairs.filt <- onecut.pairs[1:30]
textplot(seq(length(onecut.pairs.filt)), onecut.pairs.filt, words = names(onecut.pairs.filt), xlab = "Index", ylab = "Number of pairs in foreground", cex = 0.9)
