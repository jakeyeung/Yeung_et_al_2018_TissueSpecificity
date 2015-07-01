# Differential sitecount analysis: associate with gene expression
# 2015-06-30

library(hash)
library(dplyr)
library(ggplot2)
library(corrplot)
library(PMA)  # penalized CCA, install impute on Bioconductor
library(CCA)

source("scripts/functions/LoadSitecounts.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

# Functions ---------------------------------------------------------------

PlotDifferentialSitecounts <- function(N, jgene, rhyth.tiss, flat.tiss, val = "motevo.value"){
  
  motevo.sub <- subset(N, gene == jgene)
  
  m <- dcast(data = motevo.sub, formula = motif ~ tissue, value.var = "motevo.value")
  rownames(m) <- m$motif
  m$motif <- NULL
  m <- as.matrix(m)
  # normalize by sum
#   m <- sweep(m, 2, colSums(m),"/")
#   m <- m[which(rowMeans(m) > 0), ]
  
  m.rhythmic <- m[, which(colnames(m) %in% rhyth.tiss)]
  
  # m.liver <- m[, "Liver"]
  # m.other <- m[, c("Cere", "Heart", "Kidney", "Lung", "Mus")]
  m.flat <- m[, which(colnames(m) %in% flat.tiss)]
  # m.flat <- m[, c("Liver", "Kidney")]
  if (!is.null(ncol(m.rhythmic))){
    m.rhythmic.avg <- rowMeans(m.rhythmic)
  } else {
    m.rhythmic.avg <- m.rhythmic
  }
  if (!is.null(ncol(m.flat))){
    m.flat.avg <- rowMeans(m.flat)
  } else {
    m.flat.avg <- m.flat
  }
  abs.diff <- m.rhythmic.avg - m.flat.avg
  abs.diff <- abs.diff[order(abs.diff, decreasing = TRUE)]
  # (head(sort(abs.diff, decreasing = TRUE), n = 100))
  
  par(mar=c(5, 15, 4.1, 2.1))
  barplot(abs.diff[1:50], names.arg = names(abs.diff[1:50]), las = 1, horiz = TRUE)  
}

MeanMat <- function(dat, cnames){
  # Take average of cnames and return
  if (length(cnames) > 1){
    dat <- rowMeans(dat[, cnames])
  } else {
    dat <- dat[, cnames]
  }
}

cancor.svd <- function(X, Y){
  # Do canonical correlation via SVD
  # http://www.nr.com/whp/notes/CanonCorrBySVD.pdf
  
  rank.X <- rankMatrix(X)[1]
  rank.Y <- rankMatrix(Y)[1]
  min.rank <- min(rank.X, rank.Y)
  print(paste("Rank:", min.rank))
  
  # Setp 1: SVD both X and Y
  s.X <- svd(X)
  s.Y <- svd(Y)
  
  # Step 2: Form u.xTu.y, then SVD THAT
  s.uxuy <- svd(t(s.X$u[, seq(min.rank)]) %*% s.Y$u[, seq(min.rank)]) 
  
  # Step 3: define a and b, matrices of linear combinations 
  d.X <- diag(s.X$d[seq(min.rank)])
  d.Y <- diag(s.Y$d[seq(min.rank)])
  
  
  a <- s.X$v[, seq(min.rank)] %*% solve(d.X) %*% s.uxuy$u[, seq(min.rank)]
  b <- s.Y$v[, seq(min.rank)] %*% solve(d.Y) %*% s.uxuy$v[, seq(min.rank)]
  
  # Do quick check that D == S and uTu == 1
  print("Checking D == S")
  D <- t(a) %*% t(X) %*% Y %*% b
  print(D)
  print(diag(s.uxuy$d))
  
  print("Checking uTu == 1")
  uTu <- t(a) %*% t(X) %*% X %*% a
  print(uTu)
  return(list(s.X = s.X, s.Y = s.Y, a = a, b = b))
}

# Load --------------------------------------------------------------------


# N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_matrix"
# suffix <- "sitecounts.merged.matrix"
N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_mean_matrix"
suffix <- "sitecounts.dist_filtered.mean.matrix"
N <- LoadSitecountsEncodeAll(maindir = N.dir, suffix = suffix, with.ensemblid = FALSE, rename.tissues = TRUE)  # merged by gene
dat.long <- LoadArrayRnaSeq()
load(file = "Robjs/dat.rhyth.relamp.pvalmin1e-5.pvalmax0.05.relampmax.0.1.meancutoff6.Robj", verbose = T)

N <- N %>%
  group_by(gene, tissue) %>%
  mutate(motevo.value.norm = motevo.value / sum(motevo.value))

tfs <- GetTFs()

# Explore -----------------------------------------------------------------


jgene <- "Tars"

rhyth.tiss <- c("Liver")
flat.tiss <- c("Cere", "Kidney", "Heart", "Lung", "Mus")
# rhyth.tiss <- c("Liver", "Kidney")
# flat.tiss <- c("Cere", "Heart", "Lung", "Mus")

PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
PlotDifferentialSitecounts(N, jgene, rhyth.tiss, flat.tiss, val = "motevo.value")


# Do CCA penalized --------------------------------------------------------

jgene <- "Tars"
# # jgene <- "Rgs16"
# # jgene <- "Rgs16"
# # jgene <- "Rcan1"
# jgene <- "Slc44a1"
# jgene <- "Foxa2"
# jgene <- "Rorc"
# jgene <- "Pxmp4"
# jgene <- "Agpat6"
# jgene <- "Ypel2"
# jgene <- "Ube2u"
jgene <- "Insig2"
jgene <- "Asap2"
# # jgene <- "Pnp"
# # jgene <- "Pi4k2a"
jgene <- "Ddc"
# # jgene <- "Ethe1"
# jgene <- "Hnf1a"
jgene <- "Hnf4a"

dhs.tiss <- unique(N$tissue)

X.exprs <- dcast(data = subset(dat.rhyth.relamp, gene %in% tfs & tissue %in% dhs.tiss), formula = gene ~ tissue, value.var = "int.rnaseq")
rownames(X.exprs) <- X.exprs$gene
X.exprs$gene <- NULL

# get MOTIF matrix
X.motif <- dcast(data = subset(N, gene == jgene), formula = motif ~ tissue, value.var = "motevo.value.norm")
rownames(X.motif) <- X.motif$motif
X.motif$motif <- NULL

# make X.expr and X.motif have same colnames
X.motif <- X.motif[, colnames(X.exprs)]

# center stuff
X.exprs <- t(scale(t(as.matrix(X.exprs)), center = TRUE, scale = TRUE))
X.exprs <- X.exprs[complete.cases(X.exprs), ]
X.motif <- t(scale(t(as.matrix(X.motif)), center = TRUE, scale = TRUE))
X.motif <- X.motif[complete.cases(X.motif), ]

# Plot the biplot (learn how to interpret this) for covariance in X.motif
p.motif <- prcomp(X.motif, center = TRUE, scale. = FALSE)
biplot(p.motif, main = jgene, cex = c(0.5, 1.6), pch = 20)


# Do canonical correlation ------------------------------------------------

# cc.cca <- cancor(t(X.exprs), t(X.motif))
cc.cca <- cc(t(X.exprs), t(X.motif))
cc.cca <- cc(t(X.motif), t(X.exprs))
cc.cca <- cancor(t(X.exprs), t(X.motif), xcenter = FALSE, ycenter = FALSE)

perm.out <- CCA.permute(t(X.exprs), t(X.motif), typex="standard", typez="standard", standardize = FALSE)
cc.penalized <- CCA(t(X.exprs), t(X.motif), typex="standard", typez="standard", penaltyx = perm.out$bestpenaltyx, penaltyz = perm.out$bestpenaltyz, standardize = FALSE, K = 2)

rownames(cc.penalized$u) <- rownames(X.exprs)
rownames(cc.penalized$v) <- rownames(X.motif)

plot(cc.penalized$u[, 1], cc.penalized$u[, 2], main = jgene)
text(cc.penalized$u[, 1], cc.penalized$u[, 2], labels = rownames(cc.penalized$u))

plot(cc.penalized$v[, 1], cc.penalized$v[, 2], main = jgene)
text(cc.penalized$v[, 1], cc.penalized$v[, 2], labels = rownames(cc.penalized$v))


# Do canonical correlations by SVD ----------------------------------------

cancor.out <- cancor.svd(t(X.exprs), t(X.motif))
rownames(cancor.out$a) <- rownames(X.exprs)
rownames(cancor.out$b) <- rownames(X.motif)

plot(cancor.out$a[, 1], cancor.out$a[, 2])
text(cancor.out$a[, 1], cancor.out$a[, 2], labels = rownames(cancor.out$a))
abline(v = 0)
abline(h = 0)

plot(cancor.out$b[, 1], cancor.out$b[, 2])
text(cancor.out$b[, 1], cancor.out$b[, 2], labels = rownames(cancor.out$b))
abline(v = 0)
abline(h = 0)


# Canonical correlation ---------------------------------------------------

# do it by hand
CanonicalCorrelation <- function(X, Y){
  S.xy <- cov(X, Y)
  S.xx <- var(X)
  S.yx <- cov(Y, X)
  S.yy <- var(Y)
  A <- eigen(solve(S.xx) %*% S.xy %*% solve(S.yy) %*% S.yx)$vectors
  B <- eigen(solve(S.yy) %*% S.yx %*% solve(S.xx) %*% S.xy)$vectors
  R <- sqrt(eigen(solve(S.yy) %*% S.yx %*% solve(S.xx) %*% S.xy)$values)
  return(list(A = A, B = B, R = R))
}

CanonicalCorrelation(t(X.exprs), t(X.motif))

# http://www.statpower.net/Content/312/Lecture%20Slides/CanonicalCorrelation.pdf
source("http://www.statpower.net/R312/Steiger R Library Functions.txt")
source("http://www.statpower.net/R312/Data 1.txt")

cc.cca <- cc(X, Y)
# cc.cca <- cancor(t(X), t(Y), xcenter = FALSE, ycenter = FALSE)
cc.cca <- cancor(X, Y, xcenter = FALSE, ycenter = FALSE)
cc.jake <- CanonicalCorrelation(X, Y)

## signs of results are random
pop <- LifeCycleSavings[, 2:3]
oec <- LifeCycleSavings[, -(2:3)]
cancor(pop, oec)


mm <- read.csv("http://www.ats.ucla.edu/stat/data/mmreg.csv")
colnames(mm) <- c("Control", "Concept", "Motivation", "Read", "Write", "Math", 
                  "Science", "Sex")
summary(mm)

psych <- mm[, 1:3]
acad <- mm[, 4:8]
cc1 <- cc(psych, acad)
