# 2016-01-11
# Redo liver and flat stuff again 

library(hash)
library(dplyr)
library(reshape2)
library(ggplot2)
library(penalizedLDA)

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FisherTestSitecounts.R")
# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.collapse.flat.dist1000.Robj", verbose=T)
load("Robjs/S.collapse.liver.dist1000.Robj", verbose=T)
load("Robjs/N.long.flatliver_genes.all_motifs.dist1000.Robj", verbose=T)


# Find motifs enriched in Liver peaks vs Flat peaks in Liver-module -------

signif.genes <- subset(fits.best, model == "Liver" & amp.avg > 0.5)$gene
# signif.genes <- subset(fits.best, model == "Liver")$gene
N.subliv <- subset(N.long.all, peak %in% S.collapse$peak & model == "Liver" & gene %in% signif.genes)
# N.subliv <- subset(N.long.all, peak %in% S.collapse.flat$peak & model == "Flat")
N.subflat <- subset(N.long.all, peak %in% S.collapse.flat$peak & model == "Flat")

N.subliv <- N.subliv %>%
  group_by(motif, peak) %>%
  summarise(sitecount = sum(sitecount)) %>%
  mutate(peak.type = "Liver")

N.subflat <- N.subflat %>%
  group_by(motif, peak) %>%
  summarise(sitecount = sum(sitecount)) %>%
  mutate(peak.type = "Liver")


peaktype.hash <- hash(as.character(S.collapse$peak), as.character(S.collapse$peak.type))
peaktype.flat.hash <- hash(as.character(S.collapse.flat$peak), as.character(S.collapse.flat$peak.type))

N.subliv$peak.type <- sapply(as.character(N.subliv$peak), function(p) peaktype.hash[[p]])
N.subliv$model <- "Liver"
N.subflat$peak.type <- sapply(as.character(N.subflat$peak), function(p) peaktype.flat.hash[[p]])
N.subflat$model <- "Flat"


# Do Fisher ---------------------------------------------------------------

# compare Liver peaks versus Flat peaks on same genes
tests.motif.liver_module.livervflat <- N.subliv  %>%
  group_by(motif) %>%
  do(FisherTestSitecounts(., cutoff = 0.5, sitecount.col = "sitecount", model.col = "peak.type"))

# compare Liver peaks in Liver-rhythmic genes versus Liver peaks in flat genes
N.sub_liv.v.flat <- rbind(subset(N.subliv, peak.type == "Liver"), subset(N.subflat, peak.type == "Liver"))

tests.motif.liver_peaks.livervflat <- N.sub_liv.v.flat %>%
  group_by(motif) %>%
  do(FisherTestSitecounts(., cutoff = 0.5, sitecount.col = "sitecount", model.col = "model"))

ggplot(tests.motif.liver_module.livervflat, aes(x = odds.ratio, y = -log10(p.value), label = motif)) + geom_text() + theme_bw(24)
ggplot(tests.motif.liver_peaks.livervflat, aes(x = odds.ratio, y = -log10(p.value), label = motif)) + geom_text() + theme_bw(24)

qplot(x = tests.motif.liver_module.livervflat$odds.ratio, y = tests.motif.liver_peaks.livervflat$odds.ratio, label = tests.motif.liver_peaks.livervflat$motif) + 
  geom_text() + geom_vline(aes(xintercept = 1)) + geom_hline(aes(yintercept = 1)) + theme_bw(24)

FisherTestSitecounts(subset(N.subliv, motif == "TCF4_dimer.p2"), cutoff = 0.6, sitecount.col = "sitecount", model.col = "peak.type", show.table = TRUE)


# Use LDA  ----------------------------------------------------------------

X <- dcast(N.subliv, formula = peak ~ motif, value.var = "sitecount", fill = 0)
rownames(X) <- X$peak; X$peak <- NULL

# add liver peaks in flat genes
X.flat <- dcast(subset(N.subflat, peak.type == "Liver"), formula = peak ~ motif, value.var = "sitecount", fill = 0)
rownames(X.flat) <- X.flat$peak; X.flat$peak <- NULL

y <- sapply(rownames(X), function(p){
  ptype <- peaktype.hash[[p]] 
  if (ptype == "Liver"){
    return(1)
  } else if (ptype == "Flat"){
    return(2)
  }
})

# append y.flat to y and rbind X
y.flat <- rep(2, nrow(X.flat))
y <- c(y, y.flat)
X <- rbind(X, X.flat)

cv.out <- PenalizedLDA.cv(X,y,type="standard",lambdas=c(1e-4,1e-3,1e-2,.1, 0.2, 0.3))
PenalizedLDA(X, y, type = "standard", lambda = 0.1, K = 1)
print(cv.out)
plot(cv.out)

# make 80% of X to be training, the rest is test
xtraini1 <- sample(x = which(y == 1), size = 0.8 * length(which(y == 1)), replace = FALSE)
xtraini2 <- sample(x = which(y == 2), size = 0.8 * length(which(y == 2)), replace = FALSE)
xtesti1 <- which(y == 1)[!which(y == 1) %in% xtraini1]
xtesti2 <- which(y == 2)[!which(y == 2) %in% xtraini2]

X.train <- rbind(X[xtraini1, ], X[xtraini2, ])
y.train <- c(rep(1, length(xtraini1)), rep(2, length(xtraini2)))
X.test <- rbind(X[xtesti1, ], X[xtesti2, ])
y.test <- c(rep(1, length(xtesti1)), rep(2, length(xtesti2)))

out <- PenalizedLDA(x = X.train, y = y.train, xte = X.test, type = "standard", lambda = 0.1, K = 1, standardized = FALSE, trace = TRUE)

# out <- PenalizedLDA(x = X, y = y, type = "standard", lambda = 0.1, K = 0.1, standardized = FALSE, trace = TRUE)
plot(x = seq(length(out$discrim)), y = out$discrim)
text(x = seq(length(out$discrim)), y = out$discrim, labels = colnames(X))


# Add cross products ------------------------------------------------------

top.n <- 10
hits <- colnames(X)[order(out$discrim)][1:top.n]

# hits <- c("KLF4.p3")
# hits <- c("RORA.p2", "RXRG_dimer.p3")
hits <- c("RORA.p2")
hits <- c("ONECUT1,2.p2")
hits <- c("NR6A1.p2")
hits <- c("RXRG_dimer.p3")
# hits <- c("ATF4.p2")

top.n <- length(hits)

X <- as.matrix(X)
X.toappend <- matrix(data = NA, nrow = nrow(X), ncol = ncol(X) * top.n, dimnames = list(rownames(X), rep(colnames(X), top.n)))
for (i in seq(top.n)){
  mat.i.start <- (i - 1) * ncol(X) + 1
  mat.i.end <- mat.i.start + ncol(X) - 1
  hit <- hits[i]
  print(hit)
  x <- X[, hit]
  Xx <- sweep(X, MARGIN = 1, STATS = x, FUN = "*")
  X.toappend[, mat.i.start:mat.i.end] <- Xx
  # add colnames
  colnames(X.toappend)[mat.i.start:mat.i.end] <- paste(hit, colnames(X), sep = "-")
}

# rdo LDA
X.cross <- cbind(X, X.toappend)

# use only hits
X.cross <- X.cross[, grepl(paste(paste0(hits, "-"), collapse = "|"), colnames(X.cross))]

out.cross <- PenalizedLDA(X.cross, y, type = "standard", lambda = 0.05, K = 1, standardized = FALSE)

plot(x = seq(length(out.cross$discrim)), y = out.cross$discrim)
text(x = seq(length(out.cross$discrim)), y = out.cross$discrim, labels = colnames(X.cross))


# Test on Adr again  ------------------------------------------------------

load("Robjs/N.long.promoters_500.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)

adr.genes <- subset(fits.best, model == "Adr")$gene
flat.genes <- subset(fits.best, model == "")$gene
rhyth.genes <- subset(fits.best, n.rhyth >= 8)$gene

N.adr <- subset(N.long, gene %in% adr.genes)
N.flat <- subset(N.long, gene %in% flat.genes)
N.rhyth <- subset(N.long, gene %in% rhyth.genes)

# which motifs separate between Adr and Flat?
N.adr.mat <- dcast(data = N.adr, formula = gene.uniq ~ motif, value.var = "sitecount", fill = 0)
N.flat.mat <- dcast(data = N.flat, formula = gene.uniq ~ motif, value.var = "sitecount", fill = 0)
N.rhyth.mat <- dcast(data = N.rhyth, formula = gene.uniq ~ motif, value.var = "sitecount", fill = 0)

labs <- c(rep(1, nrow(N.adr.mat)), rep(2, nrow(N.flat.mat)))

N.adrflat <- rbind(N.adr.mat, N.flat.mat)

rownames(N.adrflat) <- N.adrflat$gene.uniq; N.adrflat$gene.uniq <- NULL

# remove column sums == 0
N.adrflat[[names(which(colSums(N.adrflat) == 0))]] <- NULL


out.adr <- PenalizedLDA(N.adrflat, labs, lambda=0.1, K=1)
# out <- PenalizedLDA(N.adrflat, labs, lambda=0, K=1)
print(out.adr)
plot(out.adr$discrim[order(out.adr$discrim)], xlim = c(-10, 200))
text(out.adr$discrim[order(out.adr$discrim)], labels = colnames(N.adrflat)[order(out.adr$discrim)])

colnames(N.adrflat)[order(out.adr$discrim, decreasing = TRUE)]
plot(out.adr$xproj, col = labs)
dat <- data.frame(xproj = out.adr$xproj, label = as.factor(labs))

ggplot(dat, aes(x = label, y = xproj)) + geom_boxplot()

# Toy ---------------------------------------------------------------------


set.seed(1)
n <- 20
p <- 100
x.toy <- matrix(rnorm(n*p), ncol=p)
y.toy <- c(rep(1,5),rep(2,5),rep(3,10))
x.toy[y.toy==1,1:10] <- x.toy[y.toy==1,1:10] + 2
x.toy[y.toy==2,11:20] <- x.toy[y.toy==2,11:20] - 2
out.toy <- PenalizedLDA(x.toy,y.toy,lambda=0.1,K=2)
print(out.toy)
plot(out.toy$xproj)
