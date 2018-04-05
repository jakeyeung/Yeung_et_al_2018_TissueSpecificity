# Jake Yeung
# Find clusters of genes by clustering log likelihood ratios yo

library(mclust)
library(reshape2)
library(cluster) 
library(fpc)  # tediously long install

# Functions ---------------------------------------------------------------

Likelihood <- function(sse, variance, Log = TRUE){
  # given sse, variance compute likelihood assuming Gaussian noise
  chi_sqr <- sse / variance
  L <- (1 / (2 * pi * variance)) * exp(-chi_sqr / 2)
  if (Log){
    L <- log(L)
  }
  return(L)
}

RatioTest <- function(l.full, l.null, Log = TRUE){
  # Given two likelihoods (full and null), calculate ratio test
  if (Log){
    D <- 2 * (l.full - l.null)
  } else {
    D <- 2 * log(l.full / l.null)
  }
}

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")

source("scripts/functions/GetClockGenes.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotFunctions.R")

# Calculate log likelihood ratio ------------------------------------------

clockgenes <- GetClockGenes()

dat.fit.likelihood <- FitRhythmicDatLong(dat.long, jget.residuals = TRUE) 

dat.fit.likelihood$l1 <- mapply(Likelihood, dat.fit.likelihood$ssq.residuals, dat.fit.likelihood$variance.1, MoreArgs = list(Log = TRUE))
dat.fit.likelihood$l2 <- mapply(Likelihood, dat.fit.likelihood$ssq.residuals.flat, dat.fit.likelihood$variance.1, MoreArgs = list(Log = TRUE))
dat.fit.likelihood$D <- mapply(RatioTest, dat.fit.likelihood$l1, dat.fit.likelihood$l2, MoreArgs = list(Log = TRUE))

expressed.genes <- subset(dat.fit.likelihood, int.rnaseq >= 4)$gene

# save(dat.fit.likelihood, file = "Robjs/dat.fit.likelihood")
load(file = "Robjs/dat.fit.likelihood")
# Cluster -----------------------------------------------------------------

mat <- dcast(subset(dat.fit.likelihood, gene %in% expressed.genes), gene ~ tissue, value.var = "D")
rownames(mat) <- mat$gene
mat$gene <- NULL
mat <- mat[complete.cases(mat), ]

fit <- kmeans(mat, centers = 12)
str(fit)
# fit <- Mclust(data = t(mat), G = 3:11, modelNames = "VVV")

# vary parameters for most readable graph

clusplot(mat, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions

plotcluster(mat, fit$cluster)

# Heatmap -----------------------------------------------------------------

PlotRelampHeatmap(mat, jtitle = "heatmap", blackend = 100, yellowstart = 101, maxval = 900, dist.method = "manhattan")


# Plot SVD for each cluster -----------------------------------------------

load("Robjs/dat.complex.maxexprs4.z_nr1d1_separately.Robj")  # contains Rik genes, good it matches nconds

filt.tiss <- c("WFAT")

table(fit$cluster)
custom.genes <- names(fit$cluster[which(fit$cluster == 3 | fit$cluster == 12)])
custom.genes <- names(fit$cluster[which(fit$cluster == 10)])  # also tissue-wide genes
custom.genes <- names(fit$cluster[which(fit$cluster == 11)])  # lung
custom.genes <- names(fit$cluster[which(fit$cluster == 4)])  # liver-ish
custom.genes <- names(fit$cluster[which(fit$cluster == 5)])  # liver-ish
custom.genes <- names(fit$cluster[which(fit$cluster == 6)])  # liver-ish


s.custom <- SvdOnComplex(subset(dat.complex, gene %in% custom.genes), value.var = jvalue.var)
eigens.custom <- GetEigens(s.custom, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.custom$u.plot, eigens.custom$v.plot, layout = jlayout)


