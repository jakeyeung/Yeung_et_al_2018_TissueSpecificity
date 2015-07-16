# Differential sitecount analysis: associate with gene expression
# 2015-06-30

setwd("~/projects/tissue-specificity/")

library(hash)
library(dplyr)
library(ggplot2)
library(corrplot)
library(PMA)  # penalized CCA, install impute on Bioconductor
library(reshape2)

source("scripts/functions/LoadSitecounts.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/DifferentialSitecountsFunctions.R")
source("scripts/functions/PlotFunctions.R")

# Functions ---------------------------------------------------------------



# Load --------------------------------------------------------------------


# N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_matrix"
# suffix <- "sitecounts.merged.matrix"
# N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_mean_matrix"
# suffix <- "sitecounts.dist_filtered.mean.matrix"
# N <- LoadSitecountsEncodeAll(maindir = N.dir, suffix = suffix, with.ensemblid = FALSE, rename.tissues = TRUE)  # merged by gene
# N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_mean_matrix_bugfixed_redo"
N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_sum_matrix_bugfixed"
suffix <- "dist_filt.bugfixed.sitecount.matrix"
N <- LoadSitecountsEncodeAll(maindir = N.dir, tissues = c("Liver", "Kidney", "Cere", "Lung", "Heart", "Mus"),
                             suffix = suffix, with.ensemblid = FALSE, rename.tissues = FALSE)  # merged by gene
dat.long <- LoadArrayRnaSeq()
load(file = "Robjs/dat.rhyth.relamp.pvalmin1e-5.pvalmax0.05.relampmax.0.1.meancutoff6.Robj", verbose = T)

# N <- N %>%
#   group_by(gene, tissue) %>%
#   mutate(motevo.value.norm = motevo.value / sum(motevo.value))

dhs.tiss <- unique(N$tissue)
tfs <- GetTFs()

# Explore -----------------------------------------------------------------


# jgene <- "Tars"
# 
# rhyth.tiss <- c("Liver")
# flat.tiss <- c("Cere", "Kidney", "Heart", "Lung", "Mus")
# # rhyth.tiss <- c("Liver", "Kidney")
# # flat.tiss <- c("Cere", "Heart", "Lung", "Mus")
# 
# PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
# PlotDifferentialSitecounts(N, jgene, rhyth.tiss, flat.tiss, val = "motevo.value")
# 
X.exprs <- dcast(data = subset(dat.rhyth.relamp, gene %in% tfs & tissue %in% dhs.tiss), formula = gene ~ tissue, value.var = "int.rnaseq")
PlotDiagnostics(N, X.exprs, dat.rhyth.relamp, "", jscale = FALSE)

# Do PCA and canonical correlation ----------------------------------------


gene.dir <- "/home/yeung/projects/tissue-specificity/data/gene_lists"
gene.files <- list()

tissues <- c("Liver", "Kidney", "Lung", "Heart")
fnames <- c("Liver_rhythmic_genes.pvalmax0.05pvalmin1e-05relamp0.1.txt", "Kidney_rhythmic_genes.pvalmax0.05pvalmin1e-04relamp0.1.txt", "Lung_rhythmic_genes.pvalmax0.05pvalmin1e-04relamp0.1.txt", "Heart_rhythmic_genes.pvalmax0.05pvalmin0.001relamp0.1.txt")

gene.files <- mapply(function(tiss, fname, gene.dir, gene.files) gene.files[[tiss]] <- file.path(gene.dir, fname), 
                     tissues, fnames, MoreArgs = list(gene.dir = gene.dir, gene.files = gene.files))

outdir <- "plots/tissue_specific_rhythmic_genes"
# init X.exprs, static. But X.motif is dynamic depending on gene
X.exprs <- dcast(data = subset(dat.rhyth.relamp, gene %in% tfs & tissue %in% dhs.tiss), formula = gene ~ tissue, value.var = "int.rnaseq")
rownames(X.exprs) <- X.exprs$gene
X.exprs$gene <- NULL

for (tiss in tissues){
  jfile <- gene.files[[tiss]]
  pdf(file.path(outdir, paste0(tiss, "_motif_variation_analysis.bugfixed.unscaled.pdf")))
  gene.list <- GetGeneListAndOrder(jfile, dat.rhyth.relamp)
  for (g in gene.list){
    PlotDiagnostics(N, X.exprs, dat.rhyth.relamp, g, jscale = FALSE)
  }
  dev.off()
}


# # Try some for fun --------------------------------------------------------
#
# library(CCA)
jgene <- "Tars"
# jscale <- FALSE
# jcenter <- TRUE
# X.motif <- SubLongToMat(N, jgene, jvar = "motevo.value")
# X.exprs.scaled <- ScaleRemoveInfs(X.exprs, center = jcenter, scale = jscale)
# X.motif.scaled <- ScaleRemoveInfs(X.motif, center = jcenter, scale = jscale)
# X.motif.exprs <- MatchColumns(X.motif.scaled, X.exprs.scaled)
# 
# cca.out <- cancor.svd(t(X.motif.exprs$X.motif), t(X.motif.exprs$X.exprs))
# 
# u <- t(X.motif.exprs$X.motif) %*% cca.out$a
# v <- t(X.motif.exprs$X.exprs) %*% cca.out$b
# 
# biplot(x = cca.out$a[, 1:2], y = u[, 1:2], main = jgene)
# biplot(x = cca.out$b[, 1:2], y = v[, 1:2], main = jgene)
# 
# # with R package
# cca.pkg <- cancor(t(X.motif.exprs$X.motif), t(X.motif.exprs$X.exprs))

# plot(cca.out$a[, 1], cca.out$a[, 2])
# text(cca.out$a[, 1], cca.out$a[, 2], labels = rownames(cca.out$a))
# abline(v = 0); abline(h = 0)
# 
# plot(cca.out$b[, 1], cca.out$b[, 2])
# text(cca.out$b[, 1], cca.out$b[, 2], labels = rownames(cca.out$b))
# abline(v = 0); abline(h = 0)
# 
# # Plot diagnostics
jgene <- "Tars"
# # jgene <- "Insig2"
# # jgene <- "Slc45a3"
# # jgene <- "Celsr1"
# # # jgene <- "Ddc"
# # # jgene <- "Polr2l"
# # jgene <- "Mfsd2a"
# # jgene <- "Camk1d"
# # jgene <- "Tjp3"
# # jgene <- "Gale"
# # jgene <- "Tars"
# jgene <- "S100a10"
# jgene <- "Mad2l2"
subset(N, gene == jgene & motif == "RORA.p2")
# subset(N, gene == jgene & motif == "bHLH_family.p2")
# # # subset(N, gene == jgene & motif == "ONECUT1.2.p2")
# # # subset(N, gene == jgene & motif == "ATF2.p2")
# # # subset(N, gene == jgene & motif == "HNF4A_NR2F1.2.p2")
# # # subset(N, gene == jgene & motif == "KLF4.p3")
# # # subset(N, gene == jgene & motif == "ATF5.p2")
# # # PlotDifferentialSitecounts(N, jgene, rhyth.tiss=c("Liver"), flat.tiss=c("Heart", "Kidney", "Lung", "Mus", "Cere"), val = "motevo.value")
# PlotDiagnostics(N, X.exprs, dat.rhyth.relamp, jgene, jscale = TRUE)
# PlotDiagnostics(N, X.exprs, dat.rhyth.relamp, jgene, jscale = FALSE)

# pdf("plots/tissue_specific_rhythmic_genes/liver_motif_variation_analysis.pdf")
# for (g in liver.genes){
#   print(g)
#   PlotDiagnostics(N, dat.rhyth.relamp, g, jscale = FALSE)
# }
# dev.off()

# jgene <- "Tars"
# # # jgene <- "Rgs16"
# # # jgene <- "Rgs16"
# # # jgene <- "Rcan1"
# # jgene <- "Slc44a1"
# # jgene <- "Foxa2"
# # jgene <- "Rorc"
# # jgene <- "Pxmp4"
# # jgene <- "Agpat6"
# # jgene <- "Ypel2"
# # jgene <- "Ube2u"
# jgene <- "Insig2"
# jgene <- "Asap2"
# # # jgene <- "Pnp"
# # # jgene <- "Pi4k2a"
# jgene <- "Ddc"
# # # jgene <- "Ethe1"
# # jgene <- "Hnf1a"
# jgene <- "Hnf4a"
# jgene <- "Galnt11"
# jgene <- "Scap"
# 
# 
# 
# X.exprs <- dcast(data = subset(dat.rhyth.relamp, gene %in% tfs & tissue %in% dhs.tiss), formula = gene ~ tissue, value.var = "int.rnaseq")
# rownames(X.exprs) <- X.exprs$gene
# X.exprs$gene <- NULL
# 
# # get MOTIF matrix
# X.motif <- dcast(data = subset(N, gene == jgene), formula = motif ~ tissue, value.var = "motevo.value.norm")
# rownames(X.motif) <- X.motif$motif
# X.motif$motif <- NULL
# 
# # make X.expr and X.motif have same colnames
# X.motif <- X.motif[, colnames(X.exprs)]
# 
# # center stuff
# X.exprs <- t(scale(t(as.matrix(X.exprs)), center = TRUE, scale = FALSE))
# X.exprs <- X.exprs[complete.cases(X.exprs), ]
# X.motif <- t(scale(t(as.matrix(X.motif)), center = TRUE, scale = FALSE))
# X.motif <- X.motif[complete.cases(X.motif), ]
# 
# # Plot the biplot (learn how to interpret this) for covariance in X.motif
# p.motif <- prcomp(X.motif, center = TRUE, scale. = FALSE)
# biplot(p.motif, main = jgene, cex = c(0.5, 1.6), pch = 20)
# 
# p.exprs <- prcomp(X.exprs, center = TRUE, scale. = FALSE)
# biplot(p.exprs, main = jgene, cex = c(0.5, 1.6), pch = 20)
# 
# # Do canonical correlations by SVD ----------------------------------------
# 
# cancor.out <- cancor.svd(t(X.exprs), t(X.motif))
# rownames(cancor.out$a) <- rownames(X.exprs)
# rownames(cancor.out$b) <- rownames(X.motif)
# 
# plot(cancor.out$a[, 1], cancor.out$a[, 2], pch = 20, main = jgene)
# text(cancor.out$a[, 1], cancor.out$a[, 2], labels = rownames(cancor.out$a))
# abline(v = 0)
# abline(h = 0)
# 
# plot(cancor.out$b[, 1], cancor.out$b[, 2], pch = 20, main = jgene)
# text(cancor.out$b[, 1], cancor.out$b[, 2], labels = rownames(cancor.out$b))
# abline(v = 0)
# abline(h = 0)
# 
# biplot(cancor.out$a[, c(1, 2)], cancor.out$b[, c(1, 2)], main = jgene)
# 
# 
# # Do canonical correlation ------------------------------------------------
# 
# # cc.cca <- cancor(t(X.exprs), t(X.motif))
# cc.cca <- cc(t(X.exprs), t(X.motif))
# cc.cca <- cc(t(X.motif), t(X.exprs))
# cc.cca <- cancor(t(X.exprs), t(X.motif), xcenter = FALSE, ycenter = FALSE)
# 
# perm.out <- CCA.permute(t(X.exprs), t(X.motif), typex="standard", typez="standard", standardize = FALSE)
# cc.penalized <- CCA(t(X.exprs), t(X.motif), typex="standard", typez="standard", penaltyx = perm.out$bestpenaltyx, penaltyz = perm.out$bestpenaltyz, standardize = FALSE, K = 2)
# 
# rownames(cc.penalized$u) <- rownames(X.exprs)
# rownames(cc.penalized$v) <- rownames(X.motif)
# 
# plot(cc.penalized$u[, 1], cc.penalized$u[, 2], main = jgene)
# text(cc.penalized$u[, 1], cc.penalized$u[, 2], labels = rownames(cc.penalized$u))
# 
# plot(cc.penalized$v[, 1], cc.penalized$v[, 2], main = jgene)
# text(cc.penalized$v[, 1], cc.penalized$v[, 2], labels = rownames(cc.penalized$v))
# 
# 
# # Canonical correlation ---------------------------------------------------
# 
# # do it by hand
# CanonicalCorrelation <- function(X, Y){
#   S.xy <- cov(X, Y)
#   S.xx <- var(X)
#   S.yx <- cov(Y, X)
#   S.yy <- var(Y)
#   A <- eigen(solve(S.xx) %*% S.xy %*% solve(S.yy) %*% S.yx)$vectors
#   B <- eigen(solve(S.yy) %*% S.yx %*% solve(S.xx) %*% S.xy)$vectors
#   R <- sqrt(eigen(solve(S.yy) %*% S.yx %*% solve(S.xx) %*% S.xy)$values)
#   return(list(A = A, B = B, R = R))
# }
# 
# CanonicalCorrelation(t(X.exprs), t(X.motif))
# 
# # http://www.statpower.net/Content/312/Lecture%20Slides/CanonicalCorrelation.pdf
# source("http://www.statpower.net/R312/Steiger R Library Functions.txt")
# source("http://www.statpower.net/R312/Data 1.txt")
# 
# cc.cca <- cc(X, Y)
# # cc.cca <- cancor(t(X), t(Y), xcenter = FALSE, ycenter = FALSE)
# cc.cca <- cancor(X, Y, xcenter = FALSE, ycenter = FALSE)
# cc.jake <- CanonicalCorrelation(X, Y)
# 
# ## signs of results are random
# pop <- LifeCycleSavings[, 2:3]
# oec <- LifeCycleSavings[, -(2:3)]
# cancor(pop, oec)
# 
# 
# mm <- read.csv("http://www.ats.ucla.edu/stat/data/mmreg.csv")
# colnames(mm) <- c("Control", "Concept", "Motivation", "Read", "Write", "Math", 
#                   "Science", "Sex")
# summary(mm)
# 
# psych <- mm[, 1:3]
# acad <- mm[, 4:8]
# cc1 <- cc(psych, acad)
