# 2016-09-14
# Jake Yeung
# pLDA with SQL db

rm(list=ls())

start <- Sys.time()

omega <- 2 * pi / 24
n <- 4
setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(reshape2)
library(ggplot2)
library(hash)
library(penalizedLDA)

source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")
source("/home/yeung/projects/posttranscriptional_regulation/functions/CosSineFunctions.R")
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
source("scripts/functions/ListFunctions.R")
source("scripts/functions/ModelStrToModel.R")
source("scripts/functions/SvdFunctions.R")

args <- commandArgs(trailingOnly=TRUE)
jmodels <- args[[1]]

tissues <- c("Liver", "Kidney")
rhyth.tiss <- ModelToTissue(jmodels)
flat.tiss <- tissues[which(tissues != rhyth.tiss)]
jmodels <- c(jmodels)

print(paste("Model:", jmodels))
print(paste("Rhyth tissue:", rhyth.tiss))
print(paste("Flat tissue:", flat.tiss))
if (length(rhyth.tiss) != 1 & length(flat.tiss) != 1){
  stop("Rhythmic tissue and flat tissue should be length 1")
}


distfilt <- 40000
jcutoff <- 3  # arbitrary
jcutoff.low <- 0  # arbitrary
# jcutoff <- 2  # arbitrary
# jcutoff <- 3  # arbitrary
cleanup <- FALSE
writepeaks <- TRUE
jmethod <- "g=1001"
saveplot <- TRUE
saverobj <- FALSE

if (is.na(distfilt)) stop("Distfilt must be numeric")

print(paste("Distance filter:", distfilt))
print(paste("DHS signal zscore cutoff:", jcutoff))


outdir <- "plots/penalized_lda/liver_kidney_wtko.sql"
dir.create(outdir)
outf.base <- paste0("plda.model.", jmodels, ".distfilt.", distfilt, ".cutoff.", jcutoff, ".cutofflow", jcutoff.low, ".method.", jmethod, ".pdf")
outf <- file.path(outdir, outf.base)
amp.min <- 0


outfile.robj <- paste0("Robjs/liver_kidney_atger_nestle/penalized_lda_mats.sql/penalized_lda_mats.sql.posterior.model.", jmodels[[1]], ".distfilt.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".cutoff", jcutoff, ".method.", jmethod, ".Robj")
# if (file.exists(outfile.robj)) stop(paste0(outfile.robj, " exists. Exiting"))



# Define rhythmic and tissue motifs ---------------------------------------

rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
rhyth.motifs <- c(rhyth.motifs, c("SRF"))
rhyth.motifs <- rhyth.motifs[which( ! rhyth.motifs %in% c("ATF6"))]
tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)
tissue.motifs <- c(tissue.motifs, c("ATF5_CREB3", "ATF6"))
tissue.motifs <- tissue.motifs[which( ! tissue.motifs %in% c("SRF"))]

# remove tissue motifs in rhyth
rhyth.motifs <- rhyth.motifs[which(!rhyth.motifs %in% intersect(rhyth.motifs, tissue.motifs))]

# redefine motifs based on MARA 
zscore.min <- 1.25
print(paste("Finding regulators for:", jmodels))
maraoutdir <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmodels, ".g=1001/atger_with_kidney.bugfixed")
if (!dir.exists(maraoutdir)){
  print(maraoutdir)
  stop("MARA Outdir does not exist")
}
act.s <- LoadActivitiesLongKidneyLiver(maraoutdir, shorten.motif.name = TRUE)
act.s.complex <- ProjectWithZscore(act.s, omega, n)
sig.motifs <- unique(as.character(subset(act.s.complex, zscore > zscore.min)$gene))
rhyth.motifs <- c(rhyth.motifs, sig.motifs)

# Functions ---------------------------------------------------------------

colMax <- function(dat){
  return(apply(dat, 1, max))
}

# Load --------------------------------------------------------------------

# Load up sitecounts from sql database
inf <- "/home/shared/sql_dbs/closestbed_multiple_genes.genomewide.merged.motifindexed.sqlite3"
motevo.tbl <- LoadDatabase(inf)

# from multigene_analysis.play_with_parameters.R 
if (!exists("fits.best")){
  load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
  fits.best <- fits.long.filt; rm(fits.long.filt)
  fits.best <- subset(fits.best, method == jmethod)
} 
if (!exists("dat.long")){
  load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
  dat.orig <- dat.long
  dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
  dat.long <- StaggeredTimepointsLivKid(dat.long)
} 
if (!exists("S.long")) load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
# load N.long.filt after we defined thhe genes we're looking for 



# Get genes and peaks -----------------------------------------------------


jgenes <- as.character(subset(fits.best, model %in% jmodels)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)


print(paste("Rhythmic genes:", length(jgenes)))
print(paste("Flat genes:", length(jgenes.flat)))


# get peaks near gene
S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
jpeaks.flat <- as.character(unique(S.sub.flat$peak))

print(paste("number of peaks in liver-specific rhythmic genes", length(jpeaks)))
print(paste("number of peaks in flat genes", length(jpeaks.flat)))


N.sub.lst <- expandingList()
# lapply(jgenes, function(jgene){
for (jgene in jgenes){
  N.long.filt.query <- filter(motevo.tbl, gene == jgene)  # peaks are not indexed, so dont take them
  N.sub.tmp <- collect(N.long.filt.query, n = Inf) 
  N.sub.lst$add(N.sub.tmp)
}
N.sub <- N.sub.lst$as.list()
N.sub <- bind_rows(N.sub)
print("Collected fg genes")
print(paste("Number of genes in N.sub:", length(unique(N.sub$gene))))
print(Sys.time() - start)



N.sub.flat.lst <- expandingList()
# N.sub.flat <- lapply(jgenes.flat, function(jgene.flat){
# jgenes.flat.test <- jgenes.flat[1:2000]
for (jgene.flat in jgenes.flat){
  N.long.filt.query.flat <- filter(motevo.tbl, gene == jgene.flat)  # peaks are not indexed, so dont take them
  N.sub.flat.tmp <- collect(N.long.filt.query.flat, n = Inf)
  N.sub.flat.lst$add(N.sub.flat.tmp)
}
N.sub.flat <- N.sub.flat.lst$as.list()
N.sub.flat <- bind_rows(N.sub.flat)

print("Collected bg genes")
print(paste("Number of genes in N.sub.flat:", length(unique(N.sub.flat$gene))))
print(Sys.time() - start)

# save(N.sub, N.sub.flat, file=paste0("/tmp/N.sub.", jmodels[[1]], ".queue.Robj"))
# print("Collected bg genes")
# print(Sys.time() - start)


if (saveplot){
  print(paste("Saving to", outf))
  pdf(outf)
}


# S.subs ------------------------------------------------------------------


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


# Assign peaks ------------------------------------------------------------


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


save(mat.fg, mat.bgnonliver, mat.bg, file = outfile.robj)
print(paste("Matrices saved to:", outfile.robj))

# Visualize selected peaks in heatmap -------------------------------------

# sanity check

print(ggplot(subset(S.sub, peak %in% mat.fg$peak & tissue %in% c("Liver", "Kidney")), aes(x = tissue, y = zscore)) + geom_boxplot() + ggtitle("Mat FG peaks"))
print(ggplot(subset(S.sub, peak %in% mat.bgnonliver$peak & tissue %in% c("Liver", "Kidney")), aes(x = tissue, y = zscore)) + geom_boxplot() + ggtitle("Mat BG Non Liver peaks"))
print(ggplot(subset(S.sub.flat, peak %in% mat.bg$peak & tissue %in% c("Liver", "Kidney")), aes(x = tissue, y = zscore)) + geom_boxplot() + ggtitle("Mat BG flat liver peaks"))

# Remove dupliactes in MAT.FG and MAT.BG
# why are there duplicates?
# mat.fg <- mat.fg[!duplicated(mat.fg), ]
# mat.bg <- mat.bg[!duplicated(mat.bg), ]

# Compare with 2 ----------------------------------------------------------

# against FLAT
# remove intersecting peaks
peaks.intersect <- intersect(mat.fg$peak, mat.bg$peak)
mat.fg <- subset(mat.fg, !peak %in% peaks.intersect)
mat.bg <- subset(mat.bg, !peak %in% peaks.intersect)

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

mat.fg <- subset(mat.fg, !peak %in% peaks.intersect)
mat.bg <- subset(mat.bg, !peak %in% peaks.intersect)

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



jlabs <-  c("Liv peaks rhyth", "Kid peaks rhyth", "Liv peaks flat")

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)



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

# PlotLdaOut(out.cross.rhythtiss3, jdim = 1, jtitle = "Discrim 1: single factors + tissrhyth cross", take.n = 30, from.bottom = FALSE)
# PlotLdaOut(out.cross.rhythtiss3, jdim = 2, jtitle = "Discrim 2: single factors + tissrhyth cross", take.n = 30, from.bottom = FALSE)

m3.dim1 <- SortLda(out.cross.rhythtiss3, jdim = 1)
m3.dim2 <- SortLda(out.cross.rhythtiss3, jdim = 2)

BoxplotLdaOut(out.cross.rhythtiss3, jdim = 1, horizontal = TRUE, axis.names = jlabs, jtitle = "Projection on 1st vector: liver vs kidney DHS peaks")
BoxplotLdaOut(out.cross.rhythtiss3, jdim = 2, horizontal = FALSE, axis.names = jlabs, jtitle = "Projection on 2nd vector: rhythmic vs flat DHS peaks")

vec.length <- sqrt(out.cross.rhythtiss3$discrim[, 1]^2 + out.cross.rhythtiss3$discrim[, 2]^2)
jsize.pairs <- vec.length * 5 + 0.01
plot(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], pch = ".",
     xlab = "Liver-specific DHS vs Kidney-specific DHS", ylab = "Liver-specific rhythmic DHS vs Nonrhythmic DHS",
     main="Motifs that separate between tissues (x-axis) and rhythmic in liver (y-axis)")
text(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], names(out.cross.rhythtiss3$x), cex = jsize.pairs)
abline(v = 0); abline(h = 0)

jsize.pairs.cut <- sapply(vec.length, function(jsize){
  if (jsize > 0.1){
    return(jsize)
  } else {
    return(0)
  }
})

labels <- names(out.cross.rhythtiss3$x)
labels.cut <- mapply(function(jlab, jsize){
  if (jsize <= 0){
    return("")
  } else {
    return(jlab)
  }
}, labels, jsize.pairs.cut)

dat.plot <- data.frame(x = out.cross.rhythtiss3$discrim[, 1],
                       y = out.cross.rhythtiss3$discrim[, 2],
                       motif = labels.cut,
                       vec.length = vec.length,
                       vec.length.cut = jsize.pairs.cut)
dat.labs <- subset(dat.plot, vec.length.cut > 0)

m <- ggplot(dat.plot, aes(x = x, y = y)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(data = dat.labs, aes(x = x, y = y, label = motif), size = 2.5) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() + 
  theme(aspect.ratio = 0.25, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Motif loadings separating liver and kidney DHS peaks") + ylab("Motif loadings separating rhythmic and flat DHS peaks")
print(m)

if (saveplot){
  dev.off()
}

if (saverobj){
  save(mat.fg, mat.bgnonliver, mat.bg, file = outfile.robj)
  print(paste("Matrices saved to:", outfile.robj))
}


print(Sys.time() - start)



# Print PDFs of web browser -----------------------------------------------

if (writepeaks){
  # pick top 20ish
  set.seed(0)
  top.n <- 25
  outdir.ucsc <- file.path(outdir, paste0(jmodels, "_ucsc_snapshots.test"))
  dir.create(outdir.ucsc)
  jhash <- hash(as.character(seq(3)), c("livpeaksFG.pdf", "nonlivpeaksBG.pdf", "flatliverpeaksBG.pdf"))
  i <- 1
  for (jdat in list(mat.fg, mat.bgnonliver, mat.bg)){
    outname <- jhash[[as.character(i)]]
    peaks <- as.character(jdat$peak)
    peaks.sub <- sample(peaks, min(top.n, length(peaks)))
    jbed <- lapply(peaks.sub, CoordToBed); jbed <- do.call(rbind, jbed)
    # write beds to file
    sink(file = file.path(outdir.ucsc, paste0(sub("^([^.]*).*", "\\1", outname), ".txt")))
    for (p in peaks.sub){
      cat(p)
      cat("\n")
    }
    sink()
    bedToUCSC(jbed, outpdf = file.path(outdir.ucsc, outname), leftwindow = 500, rightwindow = 500, jtempdir = "/tmp", theURL = "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgsid=217673714_AxABFQhpgLlfTBk1BK7Kd8aaUxaI&hgt.psOutput=on")
    i <- i + 1
  }
}
