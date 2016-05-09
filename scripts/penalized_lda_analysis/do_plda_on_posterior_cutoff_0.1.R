# 2016-04-21
# Jake Yeung
# pLDA was done before with cutoff at 0.5, but no real reason to do that. Redo with 0.1!

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

args <- commandArgs(trailingOnly=TRUE)
# distfilt <- as.numeric(args[1])
# jcutoff <- as.numeric(args[2])
distfilt <- 5000
jcutoff <- 1.5  # arbitrary
cleanup <- FALSE
writepeaks <- FALSE
if (is.na(distfilt)) stop("Distfilt must be numeric")

print(paste("Distance filter:", distfilt))
print(paste("DHS signal zscore cutoff:", jcutoff))

saveplot <- FALSE
outdir <- "/home/yeung/projects/tissue-specificity/bedfiles/lda_analysis/filtered_peaks_multigene"
outf <- paste0("plots/penalized_lda/", ".2D.posterior.multigene.distfilt.", distfilt, ".cutoff.", jcutoff, ".pdf")
# colnames(N.multigene) <- c("chromo", "start", "end", "motif", "sitecount", "chromo.peak", "start.peak", "end.peak", "gene", "dist")
amp.min <- 0
rhyth.tiss <- c("Liver")
# rhyth.tiss <- c("Liver", "Kidney")
outfile.robj <- paste0("Robjs/penalized_lda_mats.posterior.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".cutoff", jcutoff, ".Robj")
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
  load("Robjs/N.long.livertwflat.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
}



# Get genes and peaks -----------------------------------------------------

jmodels <- c("Liver")
# jmodels <- c("Liver", "Adr;Liver")

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

# bind rows create labels

mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bgnonliver)
# mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bg)

mat.fgbg <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)

colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)

out <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  

BoxplotLdaOut(out, jtitle = paste0("Single factor separation. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
PlotLdaOut(out, jtitle = paste0("Single factor loadings. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))

m.singles <- SortLda(out)

# Compare with flat -------------------------------------------------------

mat.fgbg.lab.lst.flat <- SetUpMatForLda(mat.fg, mat.bg)

mat.fgbg.flat <- mat.fgbg.lab.lst.flat$mat.fgbg; labels.flat <- mat.fgbg.lab.lst.flat$labels
colnames(mat.fgbg.flat) <- sapply(colnames(mat.fgbg.flat), RemoveP2Name)
colnames(mat.fgbg.flat) <- sapply(colnames(mat.fgbg.flat), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)

out.flat <- PenalizedLDA(mat.fgbg.flat, labels.flat, lambda = 0.1, K = 1, standardized = FALSE)  

BoxplotLdaOut(out.flat, jtitle = paste0("Single factor separation vs flat. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
PlotLdaOut(out.flat, jtitle = paste0("Single factor loadings vs flat. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))

m.singles.flat <- SortLda(out.flat)

# Compare with all 3 ------------------------------------------------------

# mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg[sample(x = seq(nrow(mat.bg)), size = 700), ], has.peaks = TRUE)
mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
# mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, has.peaks = FALSE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)

jlambda <- 0.1  # liv only
out.3 <- PenalizedLDA(mat.fgbg.3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

PlotLdaOut(out.3)

print(qplot(out.3$xproj[, 1], out.3$xproj[, 2], colour = as.factor(labels3), geom = "point", alpha = I(0.2)) + ggtitle("2D plot single factors"))

PlotLdaOut(out.3, jdim = 1, jtitle = "Discrim 1: single factors", take.n = 30, from.bottom = 1)
PlotLdaOut(out.3, jdim = 2, jtitle = "Discrim 2: single factors", take.n = 30, from.bottom = 1)

PlotLdaOut2D(out.3, jcex = 0.5)

m3 <- SortLda(out.3, jdim = 1)
m3 <- SortLda(out.3, jdim = 2)

# Do cross products -------------------------------------------------------
# 
# # do crosses
# mat.fgbg.cross <- CrossProduct(mat.fgbg, remove.duplicates = TRUE)
# dim(mat.fgbg.cross)
# # add single factors
# mat.fgbg.cross <- cbind(mat.fgbg, mat.fgbg.cross)
# # remove columns with 0 variance 
# mat.fgbg.cross[which(colSums(mat.fgbg.cross) == 0)] <- list(NULL)
# 
# # jlambda <- 0.029  # kidliv
# jlambda <- 0.015  # liv only
# out.cross <- PenalizedLDA(mat.fgbg.cross, labels, lambda = jlambda, K = 1, standardized = FALSE)
# m <- SortLda(out.cross)
# print(length(m))
# BoxplotLdaOut(out.cross, jtitle = "Cross product separation")
# PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
# PlotLdaOut(out.cross, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
# 
# 
# # Do cross product with all 3 ---------------------------------------------
# 
# mat.fgbg.cross.3 <- CrossProduct(mat.fgbg.3, remove.duplicates = TRUE)
# dim(mat.fgbg.cross.3)
# # add single factors
# mat.fgbg.cross.3 <- cbind(mat.fgbg.3, mat.fgbg.cross.3)
# # remove columns with 0 variance 
# mat.fgbg.cross.3[which(colSums(mat.fgbg.cross.3) == 0)] <- list(NULL)
# labels3.to2 <- labels3; labels3.to2[which(labels3.to2 == 3)] <- 2
# 
# # jlambda <- 0.029  # kidliv
# jlambda <- 0.01  # liv only
# out.cross.3 <- PenalizedLDA(mat.fgbg.cross.3, labels3.to2, lambda = jlambda, K = 1, standardized = FALSE)
# m <- SortLda(out.cross.3)
# print(length(m))
# BoxplotLdaOut(out.cross.3, jtitle = "Cross product separation")
# PlotLdaOut(out.cross.3, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
# PlotLdaOut(out.cross.3, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
# 

# Cross prod on tiss and rhyth --------------------------------------------
rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)

# remove tissue motifs in rhyth
rhyth.motifs <- rhyth.motifs[which(!rhyth.motifs %in% intersect(rhyth.motifs, tissue.motifs))]

mat1 <- subset(mat.fgbg, select = intersect(rhyth.motifs, colnames(mat.fgbg)))
mat2 <- subset(mat.fgbg, select = intersect(tissue.motifs, colnames(mat.fgbg)))
mat12 <- CrossProductTwoSets(mat1, mat2)

mat.fgbg.cross.rhythtiss <- cbind(mat.fgbg, mat12)
# remove columns with 0 variance 
mat.fgbg.cross.rhythtiss[which(colSums(mat.fgbg.cross.rhythtiss) == 0)] <- list(NULL)

# jlambda <- 0.029  # kidliv
jlambda <- 0.015  # liv only
out.cross.rhythtiss <- PenalizedLDA(mat.fgbg.cross.rhythtiss, labels, lambda = jlambda, K = 1, standardized = FALSE)
m.rhythtiss <- SortLda(out.cross.rhythtiss)
print(length(m.rhythtiss))
BoxplotLdaOut(out.cross.rhythtiss, jtitle = "Cross product tissue and rhyth only separation")
PlotLdaOut(out.cross.rhythtiss, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings. Tiss and rhyth. (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
PlotLdaOut(out.cross.rhythtiss, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings. Tiss and rhyth. (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))


# Cross prod on tiss and rhyth: use two background sets -------------------

mat.rhyth3 <- subset(mat.fgbg.3, select = intersect(rhyth.motifs, colnames(mat.fgbg)))
mat.tiss3 <- subset(mat.fgbg.3, select = intersect(tissue.motifs, colnames(mat.fgbg)))
mat.rhythtiss3 <- CrossProductTwoSets(mat.rhyth3, mat.tiss3)

mat.fgbg.cross.rhythtiss3 <- cbind(mat.fgbg.3, mat.rhythtiss3)
# remove columns with 0 variance 
mat.fgbg.cross.rhythtiss3[which(colSums(mat.fgbg.cross.rhythtiss3) == 0)] <- list(NULL)
labels3.to2 <- labels3; labels3.to2[which(labels3.to2 == 3)] <- 2

jlambda <- 0.015  # liv only
out.cross.rhythtiss3 <- PenalizedLDA(mat.fgbg.cross.rhythtiss3, labels3.to2, lambda = jlambda, K = 1, standardized = FALSE)
m.rhythtiss3 <- SortLda(out.cross.rhythtiss3)
print(length(m.rhythtiss3))
BoxplotLdaOut(out.cross.rhythtiss3, jtitle = "Cross product tissue and rhyth only separation (label 3 converted to 2)")
PlotLdaOut(out.cross.rhythtiss3, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings. Tiss and rhyth. (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
PlotLdaOut(out.cross.rhythtiss3, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings. Tiss and rhyth. (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))


# Cross prod on tiss and rhyth: do 2D -------------------------------------

jlambda <- 0.025  # liv only
out.cross.rhythtiss3 <- PenalizedLDA(mat.fgbg.cross.rhythtiss3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

print(qplot(out.cross.rhythtiss3$xproj[, 1], out.cross.rhythtiss3$xproj[, 2], colour = as.factor(labels3), geom = "point", alpha = I(0.2)) + ggtitle("2D plot single factors + tissuerhyth cross"))

PlotLdaOut(out.cross.rhythtiss3, jdim = 1, jtitle = "Discrim 1: single factors + tissrhyth cross", take.n = 30, from.bottom = TRUE)
PlotLdaOut(out.cross.rhythtiss3, jdim = 2, jtitle = "Discrim 2: single factors + tissrhyth cross", take.n = 30, from.bottom = TRUE)

m3.dim1 <- SortLda(out.cross.rhythtiss3, jdim = 1)
m3.dim2 <- SortLda(out.cross.rhythtiss3, jdim = 2)

# PlotLdaOut2D(out.cross.rhythtiss3, jcex = 0.5)
plot(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], pch = ".")
text(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], names(out.cross.rhythtiss3$x))

# boxplots on dim1
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 1)
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 2)


# Cross product clock OR Liver --------------------------------------------

# set.seed(0)
# liv.motifs <- c("HNF1A", "FOXA2", "HNF4A_NR2F1.2", "CUX2", "CEBPA.B_DDIT3", "RXRG_dimer", "ONECUT1.2")
# liv.motifs <- c("HNF1A", "FOXA2", "HNF4A_NR2F1.2", "CUX2", "CEBPA.B_DDIT3", "RXRG_dimer", "ONECUT1.2")
liv.motifs <- strsplit("HNF1A;FOXA2;HNF4A_NR2F1.2;CUX2;CEBPA.B_DDIT3;RXRG_dimer", split = ";")[[1]]  # from optimization
# liv.motifs <- strsplit("FOXA2;CUX2;CEBPA.B_DDIT3;RXRG_dimer", split = ";")[[1]]  # from optimization
# liv.motifs <- "FOXA2"

print(liv.motifs)

# liv.motifs <- sample(x = tissue.motifs, size = 6)

# Add a Liver-motif column

# FLAT
# liver.counts <- apply(subset(mat.fgbg.flat, select = liv.motifs), 1, sum)
# mat.fgbg.liver <- cbind(mat.fgbg.flat, Liver=liver.counts)
# mat.rhyth <- subset(mat.fgbg.liver, select = intersect(rhyth.motifs, colnames(mat.fgbg.flat)))
# jlabs <- labels.flat

# NONRHYTH
  liver.counts <- apply(subset(mat.fgbg, select = liv.motifs), 1, sum)
  mat.fgbg.liver <- cbind(mat.fgbg, Liver=liver.counts)
  mat.rhyth <- subset(mat.fgbg.liver, select = intersect(rhyth.motifs, colnames(mat.fgbg)))
  jlabs <- labels

mat.liverOR <- subset(mat.fgbg.liver, select = c(Liver))

mat.liverclock <- CrossProductTwoSets(mat.rhyth, mat.liverOR)
#   mat.liverclock <- CrossProductTwoSets(mat.fgbg, mat.liverOR)

mat.fgbg.cross.rhythLiverOR <- cbind(mat.fgbg.liver, mat.liverclock)
# remove columns with 0 variance 
mat.fgbg.cross.rhythLiverOR[which(colSums(mat.fgbg.cross.rhythLiverOR) == 0)] <- list(NULL)

# jlambda <- 0.029  # kidliv
jlambda <- 0.015  # liv only
out.cross.rhythLiverOR <- PenalizedLDA(mat.fgbg.cross.rhythLiverOR, jlabs, 
                                       lambda = jlambda, K = 1, standardized = FALSE)
m.rhythLiverOR <- SortLda(out.cross.rhythLiverOR)
print(length(m.rhythLiverOR))
BoxplotLdaOut(out.cross.rhythLiverOR, jtitle = "Cross product separation. liverOR and rhyth.")
PlotLdaOut(out.cross.rhythLiverOR, take.n = 50, from.bottom = TRUE, 
           jtitle = paste0("Cross product loadings. LiverOR and rhyth. (from bottom). Dist:", 
                           distfilt, "\nN FG peaks:", 
                           length(unique(mat.fg$gene)), 
                           "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
PlotLdaOut(out.cross.rhythLiverOR, 
           take.n = 50, from.bottom = FALSE, 
           jtitle = paste0("Cross product loadings. Tiss and rhyth. (from top). Dist:", 
                           distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), 
                           "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))

# liver and clock pairs in 3D ---------------------------------------------

liver.counts.3 <- apply(subset(mat.fgbg.3, select = liv.motifs), 1, sum)
mat.fgbg.liver.3 <- cbind(mat.fgbg.3, Liver=liver.counts.3)
mat.rhyth.3 <- subset(mat.fgbg.liver.3, select = intersect(rhyth.motifs, colnames(mat.fgbg.3)))
jlabs <- labels3
convert.3.to.2 <- FALSE
# make 3 into 2 to see if that gives us only the pairs
if (convert.3.to.2){
  jlabs[which(jlabs == 3)] <- 2
}

mat.liverOR <- subset(mat.fgbg.liver.3, select = c(Liver))

mat.liverclock <- CrossProductTwoSets(mat.rhyth, mat.liverOR)
#   mat.liverclock <- CrossProductTwoSets(mat.fgbg, mat.liverOR)

mat.fgbg.cross.rhythLiverOR <- cbind(mat.fgbg.liver.3, mat.liverclock)
# remove columns with 0 variance 

bad.motifs <- list()
for (jlab in unique(jlabs)){
  #   bad.motifs[[jlab]] <- which(colSums(mat.fgbg.cross.rhythLiverOR[which(jlabs == jlab), ]) == 0)
  bad.motifs[[jlab]] <- which(apply(mat.fgbg.cross.rhythLiverOR[which(jlabs == jlab), ], 2, sd) == 0)
  # bad.motifs[[jlab]] <- which(colSums(mat.fgbg.cross.rhythLiverOR[which(jlabs == jlab), ]) == 0)
}
bad.motifs <- unlist(bad.motifs)
mat.fgbg.cross.rhythLiverOR[bad.motifs] <- list(NULL)


# jlambda <- 0.029  # kidliv
jlambda <- 0.03  # liv only
if (convert.3.to.2){
  out.cross.rhythLiverOR.3 <- PenalizedLDA(mat.fgbg.cross.rhythLiverOR, jlabs, 
                                           lambda = jlambda, K = 1, standardized = FALSE)
  m.rhythLiverOR <- SortLda(out.cross.rhythLiverOR)
  print(length(m.rhythLiverOR))
  BoxplotLdaOut(out.cross.rhythLiverOR, jtitle = "Cross product separation. liverOR and rhyth.")
  PlotLdaOut(out.cross.rhythLiverOR, take.n = 50, from.bottom = TRUE, 
             jtitle = paste0("Cross product loadings. LiverOR and rhyth. (from bottom). Dist:", 
                             distfilt, "\nN FG peaks:", 
                             length(unique(mat.fg$gene)), 
                             "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
  PlotLdaOut(out.cross.rhythLiverOR, 
             take.n = 50, from.bottom = FALSE, 
             jtitle = paste0("Cross product loadings. Tiss and rhyth. (from top). Dist:", 
                             distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), 
                             "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
} else {
  out.cross.rhythLiverOR.3 <- PenalizedLDA(mat.fgbg.cross.rhythLiverOR, jlabs, 
                                           lambda = jlambda, K = 2, standardized = FALSE)
  m.rhythLiverOR <- SortLda(out.cross.rhythLiverOR.3)
  
  PlotLdaOut(out.cross.rhythLiverOR.3)
  
  print(qplot(out.cross.rhythLiverOR.3$xproj[, 1], out.cross.rhythLiverOR.3$xproj[, 2], colour = as.factor(jlabs), geom = "point", alpha = I(0.2)) + ggtitle("2D plot single factors"))
  
  PlotLdaOut(out.cross.rhythLiverOR.3, jdim = 1, jtitle = "Discrim 1: single factors", take.n = 30, from.bottom = TRUE)
  PlotLdaOut(out.cross.rhythLiverOR.3, jdim = 2, jtitle = "Discrim 2: single factors", take.n = 30, from.bottom = TRUE)
  
  PlotLdaOut2D(out.cross.rhythLiverOR.3, jcex = 0.5)
  
  plot(out.cross.rhythLiverOR.3$discrim[, 1], out.cross.rhythLiverOR.3$discrim[, 2], pch = ".")
  text(out.cross.rhythLiverOR.3$discrim[, 1], out.cross.rhythLiverOR.3$discrim[, 2], names(out.cross.rhythLiverOR.3$x))
  abline(h = 0); abline(v = 0)
  print(length(out.cross.rhythLiverOR.3)) 
  
  # show the boxplots that discriminate 1 and 2
  BoxplotLdaOut(out.cross.rhythLiverOR.3, jdim = 1)
  BoxplotLdaOut(out.cross.rhythLiverOR.3, jdim = 2)
}




# Is this set of Liver motifs that we picked ideal? -----------------------


# Optimize by 2 labels ----------------------------------------------------

# 
# liv.motifs <- c("HNF1A", "FOXA2", "HNF4A_NR2F1.2", "CUX2", "CEBPA.B_DDIT3", "RXRG_dimer", "ONECUT1.2", "RFX1..5_RFXANK_RFXAP")
# liv.motifs <- c("HNF1A", "FOXA2", "HNF4A_NR2F1.2", "CUX2", "ONECUT1.2")
# # liv.motifs <- c("FOXA2")
# # set.seed(0)
# # liv.motifs <- sample(x = tissue.motifs, size = 7)
# 
# print(liv.motifs)
# 
# loadings <- expandingList()
# crits <- expandingList()
# crits2 <- expandingList()
# 
# for (m in seq(length(liv.motifs))){
#   print(m)
#   combos <- combn(liv.motifs, m = m)
#   for (coli in seq(ncol(combos))){
#     livm <-  combos[, coli]
#     print(livm)
#     
#     liver.counts <- apply(subset(mat.fgbg, select = livm), 1, sum)
#     # liver.counts <- apply(subset(mat.fgbg, select = livm), 1, max)
#     mat.fgbg.liver <- cbind(mat.fgbg, Liver=liver.counts)  # add liver only
#     mat.rhyth <- subset(mat.fgbg.liver, select = intersect(rhyth.motifs, colnames(mat.fgbg)))
#     jlabs <- labels
#     
#     mat.liverOR <- subset(mat.fgbg.liver, select = c(Liver))
#     mat.liverclock <- CrossProductTwoSets(mat.rhyth, mat.liverOR)
#     mat.fgbg.cross.rhythLiverOR <- cbind(mat.fgbg.liver, mat.liverclock)
#     # remove columns with 0 variance 
#     mat.fgbg.cross.rhythLiverOR[which(colSums(mat.fgbg.cross.rhythLiverOR) == 0)] <- list(NULL)
# #     mat.fgbg.cross.rhythLiverOR[which(colMax(mat.fgbg.cross.rhythLiverOR) == 0)] <- list(NULL)
#     # jlambda <- 0.029  # kidliv
#     jlambda <- 0.03  # liv only
#     out.cross.rhythLiverOR <- PenalizedLDA(mat.fgbg.cross.rhythLiverOR, jlabs, 
#                                            lambda = jlambda, K = 1, standardized = FALSE)
#     liver.loading <- out.cross.rhythLiverOR$discrim[which(names(out.cross.rhythLiverOR$x) == "Liver")]
#     motifs.name <- paste(livm, collapse = ";")
#     out.df <- data.frame("motif"=motifs.name, "loading"=liver.loading)
#     loadings$add(out.df)
#     crits$add(out.cross.rhythLiverOR$crits[[1]][length(out.cross.rhythLiverOR$crits)])
#     crits2$add(out.cross.rhythLiverOR$crits[[2]][length(out.cross.rhythLiverOR$crits)])
#   }
# }
# loadings.lst <- loadings$as.list()
# crits.lst <- crits$as.list()
# crits2.lst <- crits2$as.list()
# 
# loadings.lst <- do.call(rbind, loadings.lst)
# crits.lst <- do.call(rbind, crits.lst)
# crits2.lst <- do.call(rbind, crits2.lst)
# 
# loadings.lst <- loadings.lst[order(loadings.lst$loading), ]
# loadings.lst$indx <- seq(1:nrow(loadings.lst))
# 
# ggplot(loadings.lst, aes(label = loadings.lst$motif, x = loadings.lst$indx, y = -1 * loadings.lst$loading)) + geom_text()
# 

# Optimize by 3 labels ----------------------------------------------------


source("scripts/functions/ListFunctions.R")
liv.motifs <- c("HNF1A", "FOXA2", "HNF4A_NR2F1.2", "CUX2", "CEBPA.B_DDIT3", "RXRG_dimer", "ONECUT1.2", "RFX1..5_RFXANK_RFXAP")
# liv.motifs <- c("HNF1A", "FOXA2", "HNF4A_NR2F1.2", "CUX2", "ONECUT1.2")
# liv.motifs <- c("FOXA2")
# set.seed(0)
# liv.motifs <- sample(x = tissue.motifs, size = 7)

print(liv.motifs)

loadings <- expandingList()

for (m in seq(length(liv.motifs))){
  print(m)
  combos <- combn(liv.motifs, m = m)
  for (coli in seq(ncol(combos))){
    livm <-  combos[, coli]
    print(livm)
    
    liver.counts.3 <- apply(subset(mat.fgbg.3, select = livm), 1, sum)
    mat.fgbg.liver.3 <- cbind(mat.fgbg.3, Liver=liver.counts.3)
    mat.rhyth.3 <- subset(mat.fgbg.liver.3, select = intersect(rhyth.motifs, colnames(mat.fgbg.3)))
    jlabs <- labels3
    
    mat.liverOR.3 <- subset(mat.fgbg.liver.3, select = c(Liver))
    mat.liverclock.3 <- CrossProductTwoSets(mat.rhyth.3, mat.liverOR.3)
    mat.fgbg.cross.rhythLiverOR.3 <- cbind(mat.fgbg.liver.3, mat.liverclock.3)
    # remove columns with 0 variance 
    mat.fgbg.cross.rhythLiverOR.3[which(colSums(mat.fgbg.cross.rhythLiverOR.3) == 0)] <- list(NULL)
    jlambda <- 0.015  # liv only
    out.cross.rhythLiverOR.3 <- PenalizedLDA(mat.fgbg.cross.rhythLiverOR.3, jlabs, 
                                           lambda = jlambda, K = 2, standardized = FALSE)
    liver.loading.3 <- out.cross.rhythLiverOR.3$discrim[which(names(out.cross.rhythLiverOR.3$x) == "Liver"), 1]
    motifs.name <- paste(livm, collapse = ";")
    obj1 <- out.cross.rhythLiverOR.3$crits[[1]][length(out.cross.rhythLiverOR.3$crits)]
    obj2 <- out.cross.rhythLiverOR.3$crits[[2]][length(out.cross.rhythLiverOR.3$crits)]
    out.df <- data.frame("motif"=motifs.name, "loading"=liver.loading.3, "obj1"=obj1, "obj2"=obj2)
    loadings$add(out.df)
  }
}
loadings.lst <- loadings$as.list()

loadings.lst <- do.call(rbind, loadings.lst)

loadings.lst <- loadings.lst[order(loadings.lst$loading), ]
loadings.lst$indx <- seq(1:nrow(loadings.lst))

ggplot(loadings.lst, aes(label = loadings.lst$motif, x = loadings.lst$indx, y = -1 * loadings.lst$loading)) + geom_text()
ggplot(loadings.lst, aes(label = loadings.lst$motif, x = loadings.lst$indx, y = obj1)) + geom_text()
ggplot(loadings.lst, aes(label = loadings.lst$motif, x = loadings.lst$indx, y = obj2)) + geom_text()



if (saveplot){
  dev.off()
}