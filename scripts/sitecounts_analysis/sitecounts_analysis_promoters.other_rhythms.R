# 2016-02-13
# Jake Yeung
# Quantify significance of sitecounts using linear discriminant analysis

library(ggplot2)
library(reshape2)
library(dplyr)
library(hash)
library("penalizedLDA")

dist.ref <- 500  # 500 left and right of promoter is reference
# sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
# dist <- 500
sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_50000_dist_sum_multigene/sitecounts.50000.multigene.mat"
dist <- 50000  # needs rescaling
# sitecounts.path <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_1000_dist_sum_multigene/sitecounts.1000.multigene.mat"
# dist <- 1000

# Function ----------------------------------------------------------------

source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ShannonEntropy.R")

GetEntropyOfGeneModule <- function(dat, stat = "mean"){
  # Given data of mean expression in each tissue, compute entropy
  # expect colname to be "exprs.mean" 
  # stat: mean or median
  if (stat == "mean"){
    dat.means <- dat %>%
      group_by(tissue) %>%
      summarise(avg.exprs.mean = mean(exprs.mean))    
  } else if (stat == "median"){
    dat.means <- dat %>%
      group_by(tissue) %>%
      summarise(avg.exprs.mean = median(exprs.mean))    
  }
  # calculate entropy
  entropy <- ShannonEntropy(dat.means$avg.exprs.mean, normalize=TRUE)
  # make into dataframe so we can use do()
  return(data.frame(entropy = entropy))
}


# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")
load("Robjs/N.long.promoters_500.Robj")

source("scripts/functions/PlotGeneAcrossTissues.R")


# Format dat --------------------------------------------------------------

dat.mean.rnaseq <- subset(dat.long, tissue != "WFAT" & experiment == "rnaseq" & gene %in% fits.best$gene) %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))

# assign each gene to a model
model.hash <- hash(as.character(fits.best$gene), as.character(fits.best$model))
dat.mean.rnaseq$model <- sapply(as.character(dat.mean.rnaseq$gene), function(g) model.hash[[g]])

# Adrenal as test ---------------------------------------------------------

jtiss <- "Adr"
Adr.genes <- as.character(subset(fits.best, model == jtiss)$gene)
PlotMeanExprsOfModel(dat.mean.rnaseq, Adr.genes, jtiss)

PlotMeanExprsOfModel(dat.mean.rnaseq, as.character(subset(fits.best, model == "Lung;Adr,Aorta,Mus")$gene), "Lung;Adr,Aorta,Mus")

min.n.genes <- 25
n.tiss <- 11

dat.entropy <- dat.mean.rnaseq %>%
  group_by(model) %>%
  filter(length(gene) > min.n.genes * n.tiss) %>%
  do(GetEntropyOfGeneModule(.)) %>%
  filter(!is.na(entropy))

plot(density(dat.entropy$entropy))

# show top
dat.entropy[order(dat.entropy$entropy), ]

# plot some examples
jmodel <- "Adr,Kidney,Liver,Mus"
jmodel <- "Adr,Aorta,Heart,Lung"
jmodel <- "Adr,Heart"
jmodel <- "BFAT,Lung"
PlotMeanExprsOfModel(dat.mean.rnaseq, as.character(subset(fits.best, model == jmodel)$gene), jmodel)


# Cluster -----------------------------------------------------------------

# TODO


# Do penalized PCA --------------------------------------------------------

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

# with CV
frac.train <- 0.8
frac.test <- 1 - frac.train

cv.out <- PenalizedLDA.cv(N.adrflat, labs, lambdas=c(0.001, 0.01, 0.1, 0.11, 0.12, 0.13))
print(cv.out)
plot(cv.out)
# Perform penalized LDA #
itrain <- sample(seq(nrow(N.adrflat)), size = frac.train * nrow(N.adrflat), replace = FALSE)
itest <- which(!seq(nrow(N.adrflat)) %in% itrain)
out <- PenalizedLDA(N.adrflat[itrain, ], labs[itrain], N.adrflat[itest, ],lambda=cv.out$bestlambda,K=cv.out$bestK)
print(out)
plot(out)
print(table(out$ypred[, out$K], labs[itest]))


# Do other ones -----------------------------------------------------------

fits.bfataorta <- subset(fits.best, n.rhyth > 1 & n.rhyth < 11)
fits.bfataorta <- fits.bfataorta[grep("(;|^)Aorta.*;BFAT(;|$)", fits.bfataorta$model), ]
jmodel <- unique(as.character(fits.bfataorta$model))
# jmodel <- "Adr"
jmodel.bg <- ""
lda.out2 <- PenalizedLdaLong(fits.best, N.long, jmodel, jmodel.bg, jlambda = 0.11, K = 1)
PlotLdaOut(lda.out2, unique(N.long$motif))

# # the Aorta,BFAT antiphasic module
fits.bfataorta <- subset(fits.best, n.rhyth > 1 & n.rhyth < 11)
fits.bfataorta <- fits.bfataorta[grep("(;|^)Aorta.*;BFAT(;|$)", fits.bfataorta$model), ]
jmodel <- unique(as.character(fits.bfataorta$model))
antiphasic.genes <- subset(fits.best, model %in% jmodel)$gene
print(length(antiphasic.genes))

N.antiphasic <- subset(N.long, gene %in% antiphasic.genes)
N.antiphasic.mat <- dcast(data = N.antiphasic, formula = gene.uniq ~ motif, value.var = "sitecount", fill = 0)
N.flat.mat <- dcast(data = N.flat, formula = gene.uniq ~ motif, value.var = "sitecount", fill = 0)

labs.antiphasic <- c(rep(1, nrow(N.antiphasic.mat)), rep(2, nrow(N.flat.mat)))

N.antiphasicflat <- rbind(N.antiphasic.mat, N.flat.mat)

rownames(N.antiphasicflat) <- N.antiphasicflat$gene.uniq; N.antiphasicflat$gene.uniq <- NULL
rownames(N.antiphasic.mat) <- N.antiphasic.mat$gene.uniq; N.antiphasic.mat$gene.uniq <- NULL
rownames(N.flat.mat) <- N.flat.mat$gene.uniq; N.flat.mat$gene.uniq <- NULL

# remove column sums == 0
N.antiphasicflat[[names(which(colSums(N.antiphasicflat) == 0))]] <- NULL

lda.out <- PenalizedLDA(N.antiphasicflat, labs.antiphasic, lambda = 0.11, K = 1)
print(lda.out)
plot(lda.out)
text(x = seq(ncol(N.antiphasicflat)), y = lda.out$discrim,labels = colnames(N.antiphasicflat))
