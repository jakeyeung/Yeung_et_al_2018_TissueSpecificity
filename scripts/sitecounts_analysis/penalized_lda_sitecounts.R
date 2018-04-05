# 2016-01-09
# Jake Yeung
# Use Penalized LDA to identify motifs separating between tissue-specific rhytmic genes and flat/rhythmic genes

# install.packages("penalizedLDA")
library("penalizedLDA")
library(reshape2)

# Test package ------------------------------------------------------------

set.seed(1)
n <- 20
p <- 100
x <- matrix(rnorm(n*p), ncol=p)
y <- c(rep(1,5),rep(2,5),rep(3,10))
x[y==1,1:10] <- x[y==1,1:10] + 2
x[y==2,90:100] <- x[y==2,90:100] - 2
out <- PenalizedLDA(x,y,lambda=.14,K=2)
print(out)

plot.penlda(out)


# Test out on HNF4A and HNF1A data ----------------------------------------

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

out <- PenalizedLDA(N.adrflat, labs, lambda=0.1, K=1)
# out <- PenalizedLDA(N.adrflat, labs, lambda=0, K=1)
print(out)
plot.penlda(out)
colnames(N.adrflat)[order(out$discrim, decreasing = TRUE)]

# with CV
frac.train <- 0.8
frac.test <- 1 - frac.train

cv.out <- PenalizedLDA.cv(N.adrflat, labs, lambdas=c(0.001, 0.01, 0.1, 0.11, 0.12, 0.13, 0.141, 0.14, 0.145, 0.15, 0.16, 0.17))
print(cv.out)
plot(cv.out)
# Perform penalized LDA #
out <- PenalizedLDA(x,y,xte=xte,lambda=cv.out$bestlambda,K=cv.out$bestK)
print(out)
plot(out)
print(table(out$ypred[,out$K],yte))

