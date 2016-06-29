# 2016-04-21
# Jake Yeung
# pLDA was done before with cutoff at 0.5, but no real reason to do that. Redo with 0.1!

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

args <- commandArgs(trailingOnly=TRUE)
distfilt <- as.numeric(args[1])
jcutoff <- as.numeric(args[2])
jcutoff.low <- as.numeric(args[3])
distfilt <- 40000
jcutoff <- 3  # arbitrary
jcutoff.low <- 0  # arbitrary
# jcutoff <- 2  # arbitrary
# jcutoff <- 3  # arbitrary
cleanup <- FALSE
writepeaks <- FALSE
jmethod <- "g=1001"

if (is.na(distfilt)) stop("Distfilt must be numeric")

print(paste("Distance filter:", distfilt))
print(paste("DHS signal zscore cutoff:", jcutoff))

saveplot <- FALSE
outdir <- "plots/penalized_lda/liver_kidney_wtko"
dir.create(outdir)
outf <- paste0(outdir, "2D.posterior.multigene.distfilt.liverWTKO.", distfilt, ".cutoff.", jcutoff, ".cutofflow", jcutoff.low, ".method.", jmethod, ".pdf")
amp.min <- 0

# jmodels <- c("Kidney_SV129")
jmodels <- c("Liver_SV129")
if (jmodels == "Kidney_SV129"){
  rhyth.tiss <- c("Kidney")
  flat.tiss <- c("Liver")
} else if (jmodels == "Liver_SV129"){
  rhyth.tiss <- c("Liver")
  flat.tiss <- c("Kidney")
}
outfile.robj <- paste0("Robjs/liver_kidney_atger_nestle/penalized_lda_mats.posterior.Liver.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".cutoff", jcutoff, ".method.", jmethod, ".Robj")
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
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotUCSC.R")
source("scripts/functions/LiverKidneyFunctions.R")


# Functions ---------------------------------------------------------------

colMax <- function(dat){
  return(apply(dat, 1, max))
}

# Load --------------------------------------------------------------------

# save_N_on_posterior_cutoff_0.1.R saves Robj image. Here we laod it up
# Do Penalized LDA as before  ---------------------------------------------

# from multigene_analysis.play_with_parameters.R 
if (!exists("fits.best")){
  load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.Robj", v=T)
  fits.best <- fits.long.filt; rm(fits.long.filt)
  fits.best <- subset(fits.best, method == jmethod)
} 
if (!exists("dat.long")){
  load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)
  dat.orig <- dat.long
  dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
  dat.long <- StaggeredTimepointsLivKid(dat.long)
} 
if (!exists("S.long")) load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
if (!exists("N.long.filt")){
#   load("Robjs/N.long.livertwflat.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
  # load("Robjs/liver_kidney_atger_nestle/N.long.3wtmodules.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
  load("Robjs/liver_kidney_atger_nestle/N.long.all_genes.3wtmodules.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
}



# Get genes and peaks -----------------------------------------------------


jgenes <- as.character(subset(fits.best, model %in% jmodels)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)


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
others.i <- which(jtiss %in% flat.tiss)

S.sub.liverpeaks <- S.sub %>%
  group_by(peak, gene) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
  filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff.low)

S.sub.nonliverpeaks <- S.sub %>%
  group_by(peak, gene) %>%
  # filter(min(zscore[others.i]) > jcutoff & max(zscore[tiss.i]) > jcutoff)
  filter(min(zscore[others.i]) > 0.5 * jcutoff & max(zscore[tiss.i]) < jcutoff.low)
  # filter(min(zscore[others.i]) > jcutoff)

jtiss.flat <- levels(S.sub.flat$tissue)
if (identical(jtiss, jtiss.flat) == FALSE){
  print("This shouldnt be necessary")
  tiss.i <- which(jtiss %in% rhyth.tiss)
  others.i <- which(jtiss %in% rhyth.tiss)
}

S.sub.flat.liverpeaks <- S.sub.flat %>%
  group_by(peak, gene) %>%
  filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff.low)
  # filter(min(zscore[tiss.i]) > jcutoff)

S.sub.flat.nonliverpeaks <- S.sub.flat %>%
  group_by(peak, gene) %>%
  filter(min(zscore[others.i]) > jcutoff & max(zscore[tiss.i]) < jcutoff.low)
  # filter(min(zscore[others.i]) > jcutoff)


# Mreg should contain peaks -----------------------------------------------

jpeak <- "chr1:72173321-72173821"
subset(S.long, peak == jpeak & tissue %in% c("Liver", "Kidney"))
subset(S.sub, peak == jpeak & tissue %in% c("Liver", "Kidney"))

# Do PLDA -----------------------------------------------------------------


if (saveplot){
  print(paste("Saving to", outf))
  pdf(outf)
}

# assign peaks by removing intersections between sets
liver.peaks.fg <- as.character(unique(S.sub.liverpeaks$peak))
liver.peaks.bg <- as.character(unique(S.sub.nonliverpeaks$peak))
intersect.peaks.nonliver <- intersect(liver.peaks.fg, liver.peaks.bg)

liver.peaks.fg <- liver.peaks.fg[which(!liver.peaks.fg %in% intersect.peaks.nonliver)]
nonliver.peaks.bg <- liver.peaks.bg[which(!liver.peaks.bg %in% intersect.peaks.nonliver)]

flat.peaks.liv.bg <- as.character(unique(S.sub.flat.liverpeaks$peak))
flat.peaks.nonliv.bg <- as.character(unique(S.sub.flat.nonliverpeaks$peak))
intersect.peaks.flat <- intersect(flat.peaks.liv.bg, flat.peaks.nonliv.bg)

flat.peaks.liv.bg <- flat.peaks.liv.bg[which(!flat.peaks.liv.bg %in% intersect.peaks.flat)]

print(paste("N.peaks rhyth tiss fg", length(liver.peaks.fg)))
print(paste("N.peaks flat tiss bg", length(nonliver.peaks.bg)))

# second filter for flat genes?? is this necessary??
liver.peaks.fg <- liver.peaks.fg[which(!liver.peaks.fg %in% intersect.peaks.flat)]
# filter flat genes for intersects
flat.peaks.liv.bg <- flat.peaks.liv.bg[which(!flat.peaks.liv.bg %in% intersect.peaks.flat)]
#   liver.peaks.
print(paste("N.peaks liver peaks in flat genes", length(flat.peaks.liv.bg)))

mat.fg <- dcast(subset(N.sub, peak %in% liver.peaks.fg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bgnonliver <- dcast(subset(N.sub, peak %in% nonliver.peaks.bg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bg <- dcast(subset(N.sub.flat, peak %in% flat.peaks.liv.bg), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)

# Visualize selected peaks in heatmap -------------------------------------

# sanity check

ggplot(subset(S.sub, peak %in% mat.fg$peak & tissue %in% c("Liver", "Kidney")), aes(x = tissue, y = zscore)) + geom_boxplot() + ggtitle("Mat FG peaks")
ggplot(subset(S.sub, peak %in% mat.bgnonliver$peak & tissue %in% c("Liver", "Kidney")), aes(x = tissue, y = zscore)) + geom_boxplot() + ggtitle("Mat BG Non Liver peaks")
ggplot(subset(S.sub.flat, peak %in% mat.bg$peak & tissue %in% c("Liver", "Kidney")), aes(x = tissue, y = zscore)) + geom_boxplot() + ggtitle("Mat BG flat liver peaks")


# Compare with 2 ----------------------------------------------------------

# against FLAT
mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bg, has.peaks = TRUE)
mat.fgbg <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)
jlambda <- 0.1  # liv only
out.flat <- PenalizedLDA(mat.fgbg, labels, lambda = jlambda, K = 1, standardized = FALSE)
PlotLdaOut(out.flat, jcex = 0.5, jtitle = "Liver peaks in flat vs rhythmic")

# against NONLIVER
mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bgnonliver, has.peaks = TRUE)
mat.fgbg <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)
jlambda <- 0.1  # liv only
out.nonliv <- PenalizedLDA(mat.fgbg, labels, lambda = jlambda, K = 1, standardized = FALSE)
PlotLdaOut(out.nonliv, jcex = 0.5, jtitle = "Liver peaks in nonrhyth vs rhyth")


# Compare with all 3 ------------------------------------------------------

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)

jlambda <- 0.025  # liv only
out.cross.3 <- PenalizedLDA(mat.fgbg.3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

jlabs <-  c("Liv rhyth", "Nonrhyth", "Liv flat")
BoxplotLdaOut(out.cross.3, jdim = 2, horizontal = FALSE, axis.names = jlabs, jtitle = "Yaxis: Liv rhyth vs others")
BoxplotLdaOut(out.cross.3, jdim = 1, horizontal = TRUE, axis.names = jlabs, jtitle = "Xaxis: liver vs nonliver")

jsize.singles <- sqrt(out.cross.3$discrim[, 1]^2 + out.cross.3$discrim[, 2]^2) * 5 + 0.01
plot(out.cross.3$discrim[, 1], out.cross.3$discrim[, 2], pch = ".", 
     xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others", 
     main="Motifs that separate between tissues (x-axis) and rhythmic in liver (y-axis)")
text(out.cross.3$discrim[, 1], out.cross.3$discrim[, 2], names(out.cross.3$x), cex = jsize.singles)
abline(v = 0); abline(h = 0)


# # Cross prod on tiss and rhyth --------------------------------------------

rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)
# add from previous analysis
# rhyth.motifs <- c(rhyth.motifs, "ATF6", "ATF5_CREB3")
tissue.motifs <- c(tissue.motifs, "ATF6", "ATF5_CREB3")

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

jsize.pairs <- sqrt(out.cross.rhythtiss3$discrim[, 1]^2 + out.cross.rhythtiss3$discrim[, 2]^2) * 5 + 0.01
plot(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], pch = ".",
     xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others",
     main="Motifs that separate between tissues (x-axis) and rhythmic in liver (y-axis)")
text(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], names(out.cross.rhythtiss3$x), cex = jsize.pairs)
abline(v = 0); abline(h = 0)
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 2, horizontal = FALSE, axis.names = jlabs, jtitle = "Yaxis: Liv rhyth vs others")
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 1, horizontal = TRUE, axis.names = jlabs, jtitle = "Xaxis: liver vs nonliver")

par(mfrow=c(2,2))
jlabs <-  c("Liv rhyth", "Nonrhyth", "Liv flat")
plot(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], pch = ".",
     xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others",
     main="Motifs that separate between tissues (x-axis) and rhythmic in liver (y-axis)")
text(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], names(out.cross.rhythtiss3$x), cex = jsize.pairs)
abline(v = 0); abline(h = 0)
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 2, horizontal = FALSE, axis.names = jlabs, jtitle = "Yaxis: Liv rhyth vs others")
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 1, horizontal = TRUE, axis.names = jlabs, jtitle = "Xaxis: liver vs nonliver")
par(mfrow=c(1,1))

if (saveplot){
  dev.off()
}
# Download browser peaks --------------------------------------------------

# pick top 20ish
set.seed(0)
top.n <- 25
outdir <- paste0("/home/yeung/projects/tissue-specificity/plots/ucsc_motif_screenshots/liver_kidney_plda_peaks.stringent", ".cutoff.", jcutoff, ".cutofflow.", jcutoff.low, ".distfilt.", distfilt)
dir.create(outdir)
jhash <- hash(as.character(seq(3)), c("livpeaksFG.pdf", "nonlivpeaksBG.pdf", "flatliverpeaksBG.pdf"))
i <- 1
for (jdat in list(mat.fg, mat.bgnonliver, mat.bg)){
  outname <- jhash[[as.character(i)]]
  peaks <- as.character(jdat$peak)
  peaks.sub <- sample(peaks, min(top.n, length(peaks)))
  jbed <- lapply(peaks.sub, CoordToBed); jbed <- do.call(rbind, jbed)
  # write beds to file
  sink(file = file.path(outdir, paste0(sub("^([^.]*).*", "\\1", outname), ".txt")))
  for (p in peaks.sub){
    cat(p)
    cat("\n")
  }
  sink()
  bedToUCSC(jbed, outpdf = file.path(outdir, outname), leftwindow = 500, rightwindow = 500, theURL = "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgsid=216113768_gOhEDfRc8B232JB2uV8N5dcx6oHf&hgt.psOutput=on")
  i <- i + 1
}
# 
# # How many are enriched? Plot the data ------------------------------------
# 
# labels.rename <- hash(as.character(seq(3)), jlabs)
# 
# mat.all <- rbind(mat.fg, mat.bgnonliver, mat.bg)
# mat.all$label <- sapply(as.character(labels3), function(l) labels.rename[[l]])
# mat.melt <- melt(mat.all, id.vars = c("gene", "peak", "label"), variable.name = "motif", value.name = "motif.count")
# 
# # # remove duplicates
# # mat.melt <- subset(mat.melt, !peak %in% duplicated(peak))
# 
# ggplot(subset(mat.melt, motif == "ZNF238.p2"), aes(x = motif.count, fill = label)) + geom_density(alpha = 0.25) + scale_x_log10()
# ggplot(subset(mat.melt, motif == "FOXA2.p3"), aes(x = motif.count, fill = label)) + geom_density(alpha = 0.25) + scale_x_log10()
# ggplot(subset(mat.melt, motif == "RORA.p2"), aes(x = motif.count, fill = label)) + geom_density(alpha = 0.25) + scale_x_log10()
# ggplot(subset(mat.melt, motif == "ONECUT1,2.p2"), aes(x = motif.count, fill = label)) + geom_density(alpha = 0.25) + scale_x_log10()
# ggplot(subset(mat.melt, motif == "bHLH_family.p2"), aes(x = motif.count, fill = label)) + geom_density(alpha = 0.25) + scale_x_log10()
# ggplot(subset(mat.melt, motif == "SRY.p2"), aes(x = motif.count, fill = label)) + geom_density(alpha = 0.25) + scale_x_log10()
# ggplot(subset(mat.melt, motif == "BPTF.p2"), aes(x = motif.count, fill = label)) + geom_density(alpha = 0.25) + scale_x_log10()
# ggplot(subset(mat.melt, motif == "MEF2{A,B,C,D}.p2"), aes(x = motif.count, fill = label)) + geom_density(alpha = 0.25) + scale_x_log10()
# 
# # how many above 0.5 posterior?
# mat.npeaks <- subset(mat.melt, motif.count > 0.5) %>%
#   group_by(label) %>%
#   summarise(n.peaks = length(unique(peak)))
# n.peaks.hash <- hash(mat.npeaks$label, mat.npeaks$n.peaks)
# 
# mat.count <- subset(mat.melt, motif.count > 0.5) %>%
#   group_by(label, motif) %>%
#   summarise(count = length(motif.count), 
#             jsum = sum(motif.count))
# 
# mat.count$count.norm <- apply(mat.count, 1, function(row){
#   l <- row[[1]]
#   npeaks <- n.peaks.hash[[l]]
#   count <- as.numeric(row[[3]])
#   count.norm <- count / npeaks
# })
# 
# 
# 

