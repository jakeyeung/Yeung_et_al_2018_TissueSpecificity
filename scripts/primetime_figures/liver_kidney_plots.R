# 2016-06-29
# Jake Yeung

rm(list=ls())

setwd("~/projects/tissue-specificity")
library(dplyr)
library(ggplot2)
library(hash)
library(ggrepel)

source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")

load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.Robj")
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj")
load("Robjs/liver_kidney_atger_nestle/dat.freq.Robj")

dat.orig <- dat.long
dat.orig <- StaggeredTimepointsLivKid(dat.orig)

dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)

# filter NA changes
dat.long <- subset(dat.long, !is.na(gene))

# PCA plot ----------------------------------------------------------------

outdir <-"plots/primetime_plots_liver_kidney_wtko" 
dir.create(outdir, showWarnings = FALSE)
pdf(file.path(outdir, "liver_kidney_wtko.pdf"))

genes <- as.character(unique(dat.long$gene))
M <- dcast(subset(dat.orig, gene %in% genes), formula = gene ~ tissue + geno + time, value.var = "exprs")

rownames(M) <- M$gene
M$gene <- NULL

M <- t(scale(t(M), center = TRUE, scale = FALSE))

p <- prcomp(M, scale=FALSE, center=FALSE)


library(wordcloud)
tiss <- sapply(rownames(p$rotation), function(sname) strsplit(sname, "_")[[1]][[1]], USE.NAMES = FALSE)
genos <- sapply(rownames(p$rotation), function(sname) strsplit(sname, "_")[[1]][[2]], USE.NAMES = FALSE)
cols <- sapply(tiss, function(tt){
  if (tt == "Liver"){
    return("blue")
  } else if (tt == "Kidney"){
    return("red")
  } else {
    warning("Must be Liver or Kidney")
  }
}, USE.NAMES = FALSE)
cols.ko <- sapply(genos, function(tt){
  if (tt == "SV129"){
    return("blue")
  } else if (tt == "BmalKO"){
    return("red")
  } else {
    warning("Must be SV129 or BmalKO")
  }
}, USE.NAMES = FALSE)

# scree
pc.var <- p$sdev ^ 2 / sum(p$sdev ^ 2)
plot(pc.var, type = "o", xlim = c(1, 10), ylab = "Fraction of Variance", xlab = "Principal Component")


pc.x <- 1
pc.y <- 2

# do PCs
for (pc in seq(1)){
  pc.x <- pc
  pc.y <- pc + 1
  pc.x.var <- signif(pc.var[pc.x], 2) * 100
  pc.y.var <- signif(pc.var[pc.y], 2) * 100
  wordcloud::textplot(x = p$rotation[, pc.x], y = p$rotation[, pc.y], words = rownames(p$rotation), 
                      col = cols, xlab = paste0("PC ", pc.x, " (", pc.x.var, "% of variance)"), ylab = paste0("PC ", pc.y, " (", pc.y.var, "% of variance)"),
                      xlim = c(-0.2, 0.2), ylim = c(-0.4, 0.4))
}



# Liver modules -----------------------------------------------------------


jmeth <- "g=1001"
i <- 1
genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% c("Liver_SV129"))$gene)
# add Rgs16
genes.tw <- c(genes.tw, "Rgs16")
# genes.tw <- genes.tw[which(genes.tw != "Hmga1-rs1")]

s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = i, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)


# Motifs under liver module -----------------------------------------------

load("Robjs/liver_kidney_atger_nestle/penalized_lda_mats.posterior.model.Liver_SV129.distfilt.40000.Liver.cutoff3.method.g=1001.Robj")

jlabs <-  c("Liv peaks rhyth", "Kid peaks rhyth", "Liv peaks flat")

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)

rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
rhyth.motifs <- c(rhyth.motifs, c("SRF"))
rhyth.motifs <- rhyth.motifs[which( ! rhyth.motifs %in% c("ATF6"))]
tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)
tissue.motifs <- c(tissue.motifs, c("ATF5_CREB3", "ATF6"))
tissue.motifs <- tissue.motifs[which( ! tissue.motifs %in% c("SRF"))]

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

ggplot(dat.plot, aes(x = x, y = y)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(data = dat.labs, aes(x = x, y = y, label = motif), size = 2.5) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() + 
  theme(aspect.ratio = 0.25, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Motif loadings separating liver and kidney DHS peaks") + ylab("Motif loadings separating rhythmic and flat DHS peaks")

ggplot(dat.plot, aes(x = x, y = y)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(data = dat.labs, aes(x = x, y = y, label = motif), size = 2.5) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() + 
  theme(aspect.ratio = 0.5, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Motif loadings separating liver and kidney DHS peaks") + ylab("Motif loadings separating rhythmic and flat DHS peaks")

# qplot(x = out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], label = labels.cut, size = I(0.02)) + 
#   geom_text(size = 4) +
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() + theme(aspect.ratio = 0.5, legend.position = "none", 
#                                                                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   xlab("Liver-specific DHS vs Kidney-specific DHS") + ylab("Liver-specific rhythmic DHS vs Nonrhythmic DHS")


# Examples ----------------------------------------------------------------

PlotGeneTissuesWTKO(subset(dat.orig, gene == "Mreg"))

dev.off()

# par(mfrow=c(2,2))
# plot(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], pch = ".",
#      xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others",
#      main="Motifs that separate between tissues (x-axis) and rhythmic in liver (y-axis)")
# text(out.cross.rhythtiss3$discrim[, 1], out.cross.rhythtiss3$discrim[, 2], names(out.cross.rhythtiss3$x), cex = jsize.pairs)
# abline(v = 0); abline(h = 0)
# BoxplotLdaOut(out.cross.rhythtiss3, jdim = 2, horizontal = FALSE, axis.names = jlabs, jtitle = "Yaxis: Liv rhyth vs others")
# BoxplotLdaOut(out.cross.rhythtiss3, jdim = 1, horizontal = TRUE, axis.names = jlabs, jtitle = "Xaxis: liver vs nonliver")
# par(mfrow=c(1,1))
