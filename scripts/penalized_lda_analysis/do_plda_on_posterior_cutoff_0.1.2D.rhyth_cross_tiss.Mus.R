# 2016-04-21
# Jake Yeung
# pLDA was done before with cutoff at 0.5, but no real reason to do that. Redo with 0.1!

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

args <- commandArgs(trailingOnly=TRUE)
# distfilt <- as.numeric(args[1])
# jcutoff <- as.numeric(args[2])
distfilt <- 5000
# jcutoff <- 1.5  # arbitrary
jcutoff <- 2  # arbitrary
cleanup <- FALSE
writepeaks <- FALSE
if (is.na(distfilt)) stop("Distfilt must be numeric")

print(paste("Distance filter:", distfilt))
print(paste("DHS signal zscore cutoff:", jcutoff))

saveplot <- FALSE
outdir <- "/home/yeung/projects/tissue-specificity/bedfiles/lda_analysis/filtered_peaks_multigene_Mus"
outf <- paste0("plots/penalized_lda/", "2D.posterior.multigene.distfilt.Mus.", distfilt, ".cutoff.", jcutoff, ".pdf")
# colnames(N.multigene) <- c("chromo", "start", "end", "motif", "sitecount", "chromo.peak", "start.peak", "end.peak", "gene", "dist")
amp.min <- 0
rhyth.tiss <- c("Mus")
outfile.robj <- paste0("Robjs/penalized_lda_mats.posterior.Mus.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".cutoff", jcutoff, ".Robj")
if (file.exists(outfile.robj)) stop(paste0(outfile.robj, " exists. Exiting"))

library(dplyr)
library(reshape2)
library(ggplot2)
library(hash)
library(penalizedLDA)

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/ShannonEntropy.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/PlotUCSC.R")  # CoordToBed
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/ListFunctions.R")


# Functions ---------------------------------------------------------------

colMax <- function(dat){
  return(apply(dat, 1, max))
}

# Load --------------------------------------------------------------------

# save_N_on_posterior_cutoff_0.1.R saves Robj image. Here we laod it up

# Do Penalized LDA as before  ---------------------------------------------

# from multigene_analysis.play_with_parameters.R 
if (!exists("fits.best")) load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
if (!exists("dat.long")) load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
if (!exists("S.long")) load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
if (!exists("N.long.filt")){
#   load("Robjs/N.long.livertwflat.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
  load("Robjs/N.long.Mus_flat.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
}



# Get genes and peaks -----------------------------------------------------

jmodels <- c("Mus")

jgenes <- as.character(subset(fits.best, model %in% jmodels & amp.avg > amp.min)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)

filter.for.cc.genes <- FALSE
if (filter.for.cc.genes){
  load("Robjs/wtko.dat.wtko.hog.commongenes.Robj", v=T)
  load("Robjs/wtko.fits.all.long.commongenes.Robj", v=T)
  clock.controlled.genes <- as.character(subset(fits.all.long.wtkohog, model %in% c("WT,Liver", "WT;KO;Liver", "WT;Liver", "Liver"))$gene)
  jgenes <- intersect(jgenes, clock.controlled.genes)
  
}

print(paste("Rhythmic genes:", length(jgenes)))
print(paste("Flat genes:", length(jgenes.flat)))


# print(paste("Rhythmic genes:", length(jgenes)))
# print(paste("Flat genes:", length(jgenes.flat)))

# get peaks near gene
S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
jpeaks.flat <- as.character(unique(S.sub.flat$peak))

print(paste("number of peaks in liver-specific rhythmic genes", length(jpeaks)))
print(paste("number of peaks in flat genes", length(jpeaks.flat)))

N.sub <- subset(N.long.filt, gene %in% jgenes & peak %in% jpeaks)  # should not have any missing peaks
N.sub.flat <- subset(N.long.filt, gene %in% jgenes.flat & peak %in% jpeaks.flat)

# Clean up ram ------------------------------------------------------------
if (cleanup){
  rm(S.long, N.long.filt)
}



# Identify tissue-specific peaks ------------------------------------------


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


# Do PLDA -----------------------------------------------------------------


if (saveplot){
  print(paste("Saving to", outf))
  pdf(outf)
}

liver.peaks.fg <- as.character(unique(S.sub.liverpeaks$peak))
nonliver.peaks.bg <- as.character(unique(S.sub.nonliverpeaks$peak))
liver.peaks.bg <- as.character(unique(S.sub.flat.liverpeaks$peak))
intersect.peaks <- intersect(liver.peaks.fg, liver.peaks.bg)
intersect.peaks.nonliver <- intersect(liver.peaks.fg, nonliver.peaks.bg)
if (length(intersect.peaks.nonliver) == 0){
  warning("Warning, there are intersecting peaks between Liver FG and Non Liver BG!")
}
print(paste("N.peaks rhyth tiss fg", length(liver.peaks.fg)))
print(paste("N.peaks flat tiss bg", length(nonliver.peaks.bg)))
print(paste("N.peaks rhyth tiss bg", length(liver.peaks.bg)))

liver.peaks.fg <- liver.peaks.fg[which(!liver.peaks.fg %in% intersect.peaks)]
liver.peaks.bg <- liver.peaks.bg[which(!liver.peaks.bg %in% intersect.peaks)]
#   liver.peaks.
print(paste("N.peaks fg", length(liver.peaks.fg)))
print(paste("N.peaks bg", length(liver.peaks.bg)))


mat.fg <- dcast(subset(N.sub, peak %in% liver.peaks.fg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bgnonliver <- dcast(subset(N.sub, peak %in% nonliver.peaks.bg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bg <- dcast(subset(N.sub.flat, peak %in% liver.peaks.bg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)

# Compare with all 3 ------------------------------------------------------

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)


# Cross prod on tiss and rhyth --------------------------------------------
rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)

# remove tissue motifs in rhyth
rhyth.motifs <- rhyth.motifs[which(!rhyth.motifs %in% intersect(rhyth.motifs, tissue.motifs))]

# Cross prod on tiss and rhyth: use two background sets -------------------

mat.rhyth3 <- subset(mat.fgbg.3, select = intersect(rhyth.motifs, colnames(mat.fgbg.3)))
mat.tiss3 <- subset(mat.fgbg.3, select = intersect(tissue.motifs, colnames(mat.fgbg.3)))
mat.rhythtiss3 <- CrossProductTwoSets(mat.rhyth3, mat.tiss3)

mat.fgbg.cross.rhythtiss3 <- cbind(mat.fgbg.3, mat.rhythtiss3)
# remove columns with 0 variance 
mat.fgbg.cross.rhythtiss3[which(colSums(mat.fgbg.cross.rhythtiss3) == 0)] <- list(NULL)

# Cross prod on tiss and rhyth: do 2D -------------------------------------

jlambda <- 0.025  # liv only
out.cross.rhythtiss3 <- PenalizedLDA(mat.fgbg.cross.rhythtiss3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

print(qplot(out.cross.rhythtiss3$xproj[, 1], out.cross.rhythtiss3$xproj[, 2], colour = as.factor(labels3), geom = "point", alpha = I(0.2)) + ggtitle("2D plot single factors + tissuerhyth cross"))

PlotLdaOut(out.cross.rhythtiss3, jdim = 1, jtitle = "Discrim 1: single factors + tissrhyth cross", take.n = 30, from.bottom = TRUE)
PlotLdaOut(out.cross.rhythtiss3, jdim = 2, jtitle = "Discrim 2: single factors + tissrhyth cross", take.n = 30, from.bottom = TRUE)

m3.dim1 <- SortLda(out.cross.rhythtiss3, jdim = 1)
m3.dim2 <- SortLda(out.cross.rhythtiss3, jdim = 2)

jsize <- sqrt(out.cross.rhythtiss3$discrim[, 1]^2 + out.cross.rhythtiss3$discrim[, 2]^2) * 5 + 0.01
plot(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], pch = ".", 
     xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others", 
     main="Motifs that separate between tissues (x-axis) and rhythmic in liver (y-axis)")
text(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], names(out.cross.rhythtiss3$x), cex = jsize)
abline(v = 0); abline(h = 0)

par(mfrow=c(2,2))
jlabs <-  c("Liv rhyth", "Nonrhyth", "Liv flat")
plot(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], pch = ".", 
     xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others", 
     main="Motifs that separate between tissues (x-axis) and rhythmic in liver (y-axis)")
text(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], names(out.cross.rhythtiss3$x), cex = jsize)
abline(v = 0); abline(h = 0)
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 2, horizontal = FALSE, axis.names = jlabs, jtitle = "Yaxis: Liv rhyth vs others")
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 1, horizontal = TRUE, axis.names = jlabs, jtitle = "Xaxis: liver vs nonliver")
par(mfrow=c(1,1))

# boxplots on dim1
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 1)
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 2)

if (saveplot){
  dev.off()
}