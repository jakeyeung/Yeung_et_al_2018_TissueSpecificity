# 2016-04-04
# Jake Yeung
# multigene_analysis.R
# Assign peaks to many genes: how does this change our analysis?

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

distfilt <- 10000
jcutoff <- 3  # arbitrary
if (is.na(distfilt)) stop("Distfilt must be numeric")

amp.min <- 0
rhyth.tiss <- c("Liver")
# rhyth.tiss <- c("Liver", "Kidney")

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
source("scripts/functions/PlotUCSC.R")  # CoordToBed
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/HandleMotifNames.R")


# Function ----------------------------------------------------------------


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

N.sub <- subset(N.long.filt, peak %in% jpeaks)  # should not have any missing peaks
N.sub.flat <- subset(N.long.filt, peak %in% jpeaks.flat)



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



# Name peaks --------------------------------------------------------------


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


# Merge sitecounts near peaks ---------------------------------------------

N.sub.liver.peaks.fg <- subset(N.sub, gene %in% jgenes & peak %in% liver.peaks.fg) %>%
  group_by(gene, motif) %>%
  summarise(sitecount = max(sitecount))
  summarise(sitecount = sum(sitecount))

N.sub.nonliver.peaks.bg <- subset(N.sub, gene %in% jgenes & peak %in% nonliver.peaks.bg) %>%
  group_by(gene, motif) %>%
  summarise(sitecount = max(sitecount))
  summarise(sitecount = sum(sitecount))

N.sub.liver.peaks.bg <- subset(N.sub.flat, gene %in% jgenes.flat & peak %in% liver.peaks.bg) %>%
  group_by(gene, motif) %>%
  summarise(sitecount = max(sitecount))
  summarise(sitecount = sum(sitecount))

# plot(density(subset(N.sub.liver.peaks.fg, motif=="RORA.p2")$sitecount))
# plot(density(subset(N.sub.nonliver.peaks.bg, motif=="RORA.p2")$sitecount))



# Do penalized LDA --------------------------------------------------------


mat.fg <- dcast(N.sub.liver.peaks.fg, formula = gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bgnonliver <- dcast(N.sub.nonliver.peaks.bg, formula = gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bg <- dcast(N.sub.liver.peaks.bg, formula = gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)

mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bgnonliver, has.peaks = FALSE)
# mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bg)

mat.fgbg <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)

out <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  

BoxplotLdaOut(out, jtitle = paste0("Single factor separation. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
PlotLdaOut(out, jtitle = paste0("Single factor loadings. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))

m.singles <- SortLda(out)


# Plot most separating genes to see if they are ZT20 ----------------------

# fits.sub <- subset(fits.best, gene %in% jgenes)
# phase.hash <- hash(as.character(fits.sub$gene), fits.sub$phase.avg)

# out.df <- data.frame(gene=genes.ordered.all, phase=sapply(genes.ordered.all, function(g) phase.hash[[g]]))

# ggplot(out.df, aes(x=gene, y=phase)) + geom_point()

# genes.ordered.all <- rownames(mat.fgbg)[order(out$xproj)]
# # take fg only
# genes.ordered.all <- genes.ordered.all[which(genes.ordered.all %in% mat.fg$gene)]
# head(genes.ordered.all)
# # plot top motifs
# jgene <- "Cnot1"
# for (jgene in genes.ordered.all[1:20]){
#   PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
#   vec <- subset(mat.fg, gene==jgene, select=-gene)
#   vec <- vec[order(vec, decreasing=TRUE)]
#   vec <- vec[which(vec > 0)]
#   vec <- t(as.matrix(vec))
#   plot(vec, main=jgene)
#   text(vec, labels = rownames(vec))
# }
# 
# # Compare with flat -------------------------------------------------------
# 
# mat.fgbg.lab.lst.flat <- SetUpMatForLda(mat.fg, mat.bg, has.peaks = FALSE)
# # mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bg)
# 
# mat.fgbg.flat <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels
# 
# out.flat <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  
# 
# BoxplotLdaOut(out.flat, jtitle = paste0("Single factor separation vs flat. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
# PlotLdaOut(out.flat, jtitle = paste0("Single factor loadings vs flat. Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
# 
# m.singles <- SortLda(out.flat)


# Do cross products -------------------------------------------------------

# do crosses

mat.fgbg.cross <- CrossProduct(mat.fgbg, remove.duplicates = TRUE)
dim(mat.fgbg.cross)
# add single factors
mat.fgbg.cross <- cbind(mat.fgbg, mat.fgbg.cross)
# remove columns with 0 variance 
mat.fgbg.cross[which(colSums(mat.fgbg.cross) == 0)] <- list(NULL)

# jlambda <- 0.029  # kidliv
jlambda <- 0.01  # liv only
out.cross <- PenalizedLDA(mat.fgbg.cross, labels, lambda = jlambda, K = 1, standardized = FALSE)
m <- SortLda(out.cross)
print(length(m))
BoxplotLdaOut(out.cross, jtitle = "Cross product separation")
PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE, jtitle = paste0("Cross product loadings (from bottom). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))
PlotLdaOut(out.cross, take.n = 50, from.bottom = FALSE, jtitle = paste0("Cross product loadings (from top). Dist:", distfilt, "\nN FG peaks:", length(unique(mat.fg$gene)), "\nN BG peaks:", length(unique(mat.bgnonliver$gene))))


# Do cross products on rhythmic and tissue regulators ---------------------

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
BoxplotLdaOut(out.cross.rhythtiss, jtitle = "Cross product separation. Tiss and rhyth.")
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
