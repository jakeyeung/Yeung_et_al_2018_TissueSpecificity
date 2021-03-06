# 2016-04-17
# Jake Yeung
# multigene_analysis..play_with_parameters.R
# Assign peaks to many genes: play with parametesr such as taking MEAN rather than MAX for tissue-specific peaks
# and expanding distance filter

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

args <- commandArgs(trailingOnly=TRUE)
# distfilt <- as.numeric(args[1])
# jcutoff <- as.numeric(args[2])
distfilt <- 5000
jcutoff <- 1.5  # arbitrary
cleanup <- TRUE
writepeaks <- FALSE
if (is.na(distfilt)) stop("Distfilt must be numeric")

outdir <- "/home/yeung/projects/tissue-specificity/bedfiles/lda_analysis/filtered_peaks_multigene"
outf <- paste0("plots/penalized_lda/", ".2D.multigene.distfilt.", distfilt, ".cutoff.", jcutoff, ".pdf")
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
  

# Do penalized LDA --------------------------------------------------------


# print(paste("Saving to", outf))
# pdf(outf)

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

PlotLdaOut(out.3, jdim = 1, jtitle = "Discrim 1: single factors")
PlotLdaOut(out.3, jdim = 2, jtitle = "Discrim 2: single factors")

PlotLdaOut2D(out.3, jcex = 0.5)

m3 <- SortLda(out.3, jdim = 1)
m3 <- SortLda(out.3, jdim = 2)

  
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

  # Cross product clock OR Liver
  
  liv.motifs <- c("HNF1A", "FOXA2", "HNF4A_NR2F1.2", "CUX2", "CEBPA.B_DDIT3", "RXRG_dimer", "ONECUT1.2")
  
  # Add a Liver-motif column

  # FLAT
  liver.counts <- apply(subset(mat.fgbg.flat, select = liv.motifs), 1, sum)
  mat.fgbg.liver <- cbind(mat.fgbg.flat, Liver=liver.counts)
  mat.rhyth <- subset(mat.fgbg.liver, select = intersect(rhyth.motifs, colnames(mat.fgbg.flat)))
  jlabs <- labels.flat

  # NONRHYTH
#   liver.counts <- apply(subset(mat.fgbg, select = liv.motifs), 1, sum)
#   mat.fgbg.liver <- cbind(mat.fgbg, Liver=liver.counts)
#   mat.rhyth <- subset(mat.fgbg.liver, select = intersect(rhyth.motifs, colnames(mat.fgbg)))
#   jlabs <- labels
  
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


# Now do this in 3D! ------------------------------------------------------

liver.counts.3 <- apply(subset(mat.fgbg.3, select = liv.motifs), 1, sum)
mat.fgbg.liver.3 <- cbind(mat.fgbg.3, Liver=liver.counts.3)
mat.rhyth.3 <- subset(mat.fgbg.liver.3, select = intersect(rhyth.motifs, colnames(mat.fgbg.3)))
jlabs <- labels3

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
jlambda <- 0.015  # liv only
out.cross.rhythLiverOR.3 <- PenalizedLDA(mat.fgbg.cross.rhythLiverOR, jlabs, 
                                       lambda = jlambda, K = 2, standardized = FALSE)

PlotLdaOut(out.cross.rhythLiverOR.3, jdim = 1, jtitle = "Discrim1: includes Liver motif")
PlotLdaOut(out.cross.rhythLiverOR.3, jdim = 2, jtitle = "Discrim2: includes Liver motif")

PlotLdaOut2D(out.cross.rhythLiverOR.3, jcex = 0.5, jtitle = "2D discrim plot: includes Liver motif")

m3 <- SortLda(out.3, jdim = 1)
m3 <- SortLda(out.3, jdim = 2)


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


# dev.off()


# # Downstream analysis -----------------------------------------------------

df.sub <- melt(mat.fg, variable.name = "motif", id.vars = c("gene", "peak"), value.name = "sitecounts")
df.sub <- df.sub %>%
  group_by(gene) %>%
  filter(sitecounts > 0) %>%
  arrange(sitecounts) %>%
  mutate(indx = order(sitecounts, decreasing = TRUE))
df.sub$motif <- as.factor(sapply(as.character(df.sub$motif), function(m){
  m <- RemoveP2Name(m)
  m <- RemoveCommasBraces(m) 
}))

rora.hits <- subset(mat.fg, RORA.p2 > 0, select = c(peak, gene, RORA.p2))
liv.motifs <- c("HNF1A", "FOXA2", "HNF4A_NR2F1.2", "CUX2", "CEBPA.B_DDIT3", "RXRG_dimer", "ONECUT1.2")

jgene <- "1110012L19Rik"
jpeak <- as.character(sample(x = rora.hits$peak, size = 1))
jgene <- as.character(subset(rora.hits, peak == jpeak)$gene)
ggplot(subset(df.sub, gene == jgene & peak == jpeak & motif %in% c(liv.motifs, "RORA")), aes(x = indx, y = sitecounts, label = motif)) + geom_text() + 
  ggtitle(paste0(jgene, "\n", jpeak))
ggplot(subset(df.sub, gene == jgene & motif %in% c(liv.motifs, "RORA")), aes(x = indx, y = sitecounts, label = motif)) + geom_text() + 
  ggtitle(paste0(jgene, "\npeaks near gene"))
PlotGeneAcrossTissues(subset(dat.long, gene == jgene))


# # How many ROR motifs have a Liver-motif adjacent? How many are alone --------
# 
# df.sub.liv <- subset(df.sub, peak %in% rora.hits$peak & motif %in% liv.motifs) %>%
#   group_by(peak, gene) %>%
#   summarise(liv.sitecounts = sum(sitecounts))
# peaks.ror_with_liv <- df.sub.liv$peak
# n.peaks <- length(df.sub.liv$peak)
# total.peaks <- length(rora.hits$peak)
# 
# 
# # And the others peaks do they still contain a Liver motif? ---------------
# 
# # other motifs
# peaks.ror_no_liv <- setdiff(rora.hits$peak, peaks.ror_with_liv)
# genes.ror_no_liv <- unique(as.character(subset(df.sub, peak %in% peaks.ror_no_liv)$gene))
# 
# df.sub.liv.gene <- subset(df.sub, gene %in% genes.ror_no_liv & motif %in% liv.motifs) %>%
#   group_by(gene) %>%
#   summarise(liv.sitecounts = sum(sitecounts))
# 
# # how many genes made it?
# genes.ror_no_liv.but_distal <- as.character(df.sub.liv.gene$gene)
# 
# 
# print(paste("Number of genes with ROR but no adjacent liv motif:", length(genes.ror_no_liv)))
# print(paste("Number of genes with ROR but no adjacent liv motif but have distal liv motif (dist:):", 
#             length(genes.ror_no_liv.but_distal)))


# Finally these other 12 are due to my distance threshold being insufficient? --------
# 
# genes.ror_no_liv.no_distal <- setdiff(genes.ror_no_liv, genes.ror_no_liv.but_distal)
# 
# subset(df.sub, gene %in% genes.ror_no_liv.no_distal & motif == "RORA")
# 
# # example: Aagab distal ROR may be further than what we thought
# motifs.with.P <- unique(as.character(N.sub.gene$motif))
# liv.motifs.with.P <- motifs.with.P[grepl(paste(liv.motifs, collapse = "|"), motifs.with.P)]
# jgene <- "Aagab"
# jgene <- "Ube2u"
# dist.filt.new <- 10000
# S.sub.gene <- subset(S.long, gene == jgene & dist < dist.filt.new)
# 
# S.sub.liverpeaks.gene <- S.sub.gene %>%
#   group_by(peak) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
#   filter(min(zscore[tiss.i]) > jcutoff & mean(zscore[others.i]) < jcutoff)
# 
# peaks.gene <- as.character(unique(S.sub.liverpeaks.gene$peak))
# N.sub.gene <- subset(N.long.filt, peak %in% peaks.gene & gene == jgene)
# df.sub.gene <- subset(N.sub.gene, motif %in% liv.motifs.with.P) %>%
#   group_by(gene) %>%
#   summarise(liv.sitecounts = sum(sitecount))
# 
# # example: Ube2u has ROR motif, further away FOXA2 and HNF4A?
# jpeak <- "chr4:100479428-100479928"
# jgene <- "Ube2u"
# subset(df.sub, gene == "Ube2u" & motif %in% liv.motifs)  # signal may be obscured by Kidney...
# 
# 



# Write peaks to file -----------------------------------------------------
if (writepeaks){
  outname <- paste0("liver.peaks.fg.rhyth.tiss.2D.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".bed")
  PeaksToBedMotevo(liver.peaks.fg, file.path(outdir, outname))
  outname <- paste0("liver.peaks.bg.flat.2D.", distfilt, ".bed")
  PeaksToBedMotevo(liver.peaks.bg, file.path(outdir, outname))
  outname <- paste0("nonliver.peaks.bg.rhyth.tiss.2D.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".bed")
  
  
  # Save object -------------------------------------------------------------
  
  save(mat.fg, mat.bgnonliver, mat.bg, mat.fgbg.cross, liver.peaks.fg, nonliver.peaks.bg, liver.peaks.bg, file = outfile.robj)
}
