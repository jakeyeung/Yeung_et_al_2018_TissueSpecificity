# 2016-04-04
# Jake Yeung
# multigene_analysis.R
# Assign peaks to many genes: how does this change our analysis?

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

# args <- commandArgs(trailingOnly=TRUE)
# distfilt <- as.numeric(args[1])
distfilt <- 5000
jcutoff <- 1.5  # arbitrary
cleanup <- FALSE
writepeaks <- FALSE
if (is.na(distfilt)) stop("Distfilt must be numeric")

outdir <- "/home/yeung/projects/tissue-specificity/bedfiles/lda_analysis/filtered_peaks_multigene"
outf <- paste0("plots/penalized_lda/", "multigene.distfilt.", distfilt, ".cutoff.", jcutoff, ".pdf")
# colnames(N.multigene) <- c("chromo", "start", "end", "motif", "sitecount", "chromo.peak", "start.peak", "end.peak", "gene", "dist")
amp.min <- 0
rhyth.tiss <- c("Liver")
# rhyth.tiss <- c("Liver", "Kidney")
outfile.robj <- paste0("Robjs/penalized_lda_mats.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".cutoff", jcutoff, ".Robj")
if (file.exists(outfile.robj)) stop(paste0(outfile.robj, " exists. Exiting"))

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
source("scripts/functions/PlotUCSC.R")  # CoordToBed
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/HandleMotifNames.R")

# Function ----------------------------------------------------------------


PeaksToBedMotevo <- function(peaks, outf){
  bed <- lapply(peaks, CoordToBed)
  bed <- do.call(rbind, bed)
  bed$blank <- "."
  bed$id <- paste0("mm10_", peaks)
  outname <- paste0("rhyth.tiss.", paste(rhyth.tiss, sep = "_"), ".bed")
  write.table(bed, file = outf, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}



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


# Load --------------------------------------------------------------------


load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/N.long.multigene.distfilt.50000.Robj", v=T)
# subset(N.long.filt, 
load("Robjs/S.long.multigene.filt.50000.Robj", v=T)

# Main --------------------------------------------------------------------



if (identical(rhyth.tiss, c("Liver"))){
  jmodels <- c("Liver")
} else if (identical(rhyth.tiss, c("Liver", "Kidney"))){
  jmodels <- c("Kidney;Liver", "Kidney,Liver")  
  kidliv.genes <- fits.best$gene[sapply(fits.best$param.list, FilterKidneyLiverGenes, amp.min = 0.15)]
  kidliv.genes <- subset(fits.best, gene %in% kidliv.genes & n.rhyth <= 7)$gene
  jmodels <- as.character(subset(fits.best, gene %in% kidliv.genes)$model)
} else {
  print("Only liver or liver-kid so far")
}
jgenes <- as.character(subset(fits.best, model %in% jmodels & amp.avg > amp.min)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)

print(paste("Rhythmic genes:", length(jgenes)))
print(paste("Flat genes:", length(jgenes.flat)))

# get peaks near gene
S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
jpeaks.flat <- as.character(unique(S.sub.flat$peak))

# ggplot(S.sub, aes(x = zscore)) + geom_histogram(bins = 100) + facet_wrap(~tissue) + theme_bw(24) + geom_vline(aes(xintercept = jcutoff))
# ggplot(S.sub.flat, aes(x = zscore)) + geom_histogram(bins = 100) + facet_wrap(~tissue) + theme_bw(24) + geom_vline(aes(xintercept = jcutoff))


# Subset sitecounts (distfilt from S) -------------------------------------

N.sub <- subset(N.long.filt, gene %in% jgenes & peak %in% jpeaks)  # should not have any missing peaks
N.sub.flat <- subset(N.long.filt, gene %in% jgenes.flat & peak %in% jpeaks.flat)
# N.sub <- subset(N.long.filt, peak %in% jpeaks)  # should not have any missing peaks
# N.sub.flat <- subset(N.long.filt, peak %in% jpeaks.flat)

# rename gene names
# genes.hash <- hash(as.character(S.sub$peak), as.character(S.sub$gene))
# genes.flat.hash <- hash(as.character(S.sub.flat$peak), as.character(S.sub.flat$gene))

# N.sub$gene <- sapply(as.character(N.sub$peak), function(p) genes.hash[[p]])
# N.sub.flat$gene <- sapply(as.character(N.sub.flat$peak), function(p) genes.flat.hash[[p]])

# Clean up ram ------------------------------------------------------------
if (cleanup){
  rm(S.long, N.long.filt)
}

# Take only peaks that are "tissue-specific?" -----------------------------

# take peaks with Liver signal greater than cutoff
jtiss <- levels(S.sub$tissue)
tiss.i <- which(jtiss %in% rhyth.tiss)
others.i <- which(!jtiss %in% rhyth.tiss)

S.sub.liverpeaks <- S.sub %>%
  group_by(peak) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
  filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff)

S.sub.nonliverpeaks <- S.sub %>%
  group_by(peak) %>%
  filter(max(zscore[tiss.i]) < jcutoff & max(zscore[others.i]) > jcutoff)  # any tissue needs to be above cutoff

jtiss.flat <- levels(S.sub.flat$tissue)
if (identical(jtiss, jtiss.flat) == FALSE){
  print("This shouldnt be necessary")
  tiss.i <- which(jtiss %in% rhyth.tiss)
  others.i <- which(jtiss %in% rhyth.tiss)
}

S.sub.flat.liverpeaks <- S.sub.flat %>%
  group_by(peak) %>%
  filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff)
  

# Do penalized LDA --------------------------------------------------------


print(paste("Saving to", outf))
pdf(outf)

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
  # mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bg)
  
  mat.fgbg <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels
  colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)
  
  out <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  
  
  BoxplotLdaOut(out, jtitle = paste0("Single factor separation. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
  PlotLdaOut(out, jtitle = paste0("Single factor loadings. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
  
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
  PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))
  PlotLdaOut(out.cross, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$peak)), "\nN BG peaks:", length(unique(mat.bgnonliver$peak))))

  # Cross prod on tiss and rhyth --------------------------------------------
  rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
  tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)
  
  # remove tissue motifs in rhyth
  rhyth.motifs <- rhyth.motifs[which(!rhyth.motifs %in% intersect(rhyth.motifs, tissue.motifs))]

  cnames <- sapply(colnames(mat.fgbg), function(cname){
    return(RemoveCommasBraces(cname))
  }, USE.NAMES = FALSE)
  colnames(mat.fgbg) <- cnames
  
  mat1 <- subset(mat.fgbg, select = intersect(rhyth.motifs, colnames(mat.fgbg)))
  mat2 <- subset(mat.fgbg, select = intersect(tissue.motifs, colnames(mat.fgbg)))
  mat12 <- CrossProductTwoSets(mat1, mat2)
  
  mat.fgbg.cross.rhythtiss <- cbind(mat.fgbg, mat12)
  # remove columns with 0 variance 
  mat.fgbg.cross.rhythtiss[which(colSums(mat.fgbg.cross.rhythtiss) == 0)] <- list(NULL)
  
  # jlambda <- 0.029  # kidliv
  jlambda <- 0.015  # liv only
  out.cross.rhythtiss <- PenalizedLDA(mat.fgbg.cross.rhythtiss, labels, lambda = jlambda, K = 1, standardized = FALSE)
  m <- SortLda(out.cross.rhythtiss)
  print(length(m))
  BoxplotLdaOut(out.cross.rhythtiss, jtitle = "Cross product tissue and rhyth only separation")
  PlotLdaOut(out.cross.rhythtiss, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings. Tiss and rhyth. (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
  PlotLdaOut(out.cross.rhythtiss, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings. Tiss and rhyth. (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))

  
  # How many RORA-partners happen? ------------------------------------------
  
  n.hits <- apply(mat.fgbg.cross[seq(nrow(mat.fg)), grepl(";", colnames(mat.fgbg.cross))], 2, function(jcol) length(which(jcol > 0)))
  n.hits.single <- apply(mat.fgbg[seq(nrow(mat.fg)), ], 2, function(jcol) length(which(jcol > 0)))
  n.hits.single.bg <- apply(mat.fgbg[seq(nrow(mat.bgnonliver)) + nrow(mat.fg), ], 2, function(jcol) length(which(jcol > 0)))
  n.hits.single.all <- apply(mat.fgbg, 2, function(jcol) length(which(jcol > 0)))
  # print(head(sort(n.hits, decreasing = TRUE), n = 20))
  
  # top hits with RORA
  jmotif <- "RORA"
  rora.pairs <- sort(n.hits[grepl(jmotif, names(n.hits))], decreasing=TRUE)
  rora.pairs.filt <- rora.pairs[1:30]
  textplot(seq(length(rora.pairs.filt)), rora.pairs.filt, words = names(rora.pairs.filt), xlab = "Index", ylab = "Number of pairs in foreground", cex = 0.9, 
           main = paste0("Single RORA hits:", n.hits.single[[jmotif]], "/", nrow(mat.fg), "\nSingle RORA hits bg:", n.hits.single.bg[[jmotif]], "/", nrow(mat.bgnonliver), "\nTotal RORA hits:", n.hits.single.all[[jmotif]]))
  # 
  # onecut.pairs <- sort(n.hits[grepl("ONECUT", names(n.hits))], decreasing=TRUE)
  # onecut.pairs.filt <- onecut.pairs[1:30]
  # textplot(seq(length(onecut.pairs.filt)), onecut.pairs.filt, words = names(onecut.pairs.filt), xlab = "Index", ylab = "Number of pairs in foreground", cex = 0.9)
  
  # where are my RORA hits?
  #     mat.rora <- subset(mat.fg, select=c(peak, gene, RORA.p2)); mat.rora <- mat.rora[order(mat.rora$RORA.p2, decreasing=TRUE), ]
  #   mat.rora[1:50, ]


dev.off()

# Write peaks to file -----------------------------------------------------
if (writepeaks){
  outname <- paste0("liver.peaks.fg.rhyth.tiss.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".bed")
  PeaksToBedMotevo(liver.peaks.fg, file.path(outdir, outname))
  outname <- paste0("liver.peaks.bg.flat.", distfilt, ".bed")
  PeaksToBedMotevo(liver.peaks.bg, file.path(outdir, outname))
  outname <- paste0("nonliver.peaks.bg.rhyth.tiss.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".bed")
  
  
  # Save object -------------------------------------------------------------
  
  save(mat.fg, mat.bgnonliver, mat.bg, mat.fgbg.cross, liver.peaks.fg, nonliver.peaks.bg, liver.peaks.bg, file = outfile.robj)
}
