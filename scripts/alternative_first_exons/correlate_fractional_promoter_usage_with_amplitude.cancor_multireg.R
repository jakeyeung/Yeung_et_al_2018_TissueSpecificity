# Jake Yeung
# 2015-10-29
# Normalize amplitude by Nr1d1 DO NOT adjust for noisy genes

source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotUCSC.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")

# Functions ---------------------------------------------------------------


# Main --------------------------------------------------------------------



# load("Robjs/tpm.mr.Robj")
load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/tpm.afe.avg.Robj")
# load("Robjs/tpm.mr.weighted.Robj")
# load("Robjs/tpm.mr.weighted.shannon.Robj")
# load("Robjs/tpm.mr.weighted.Robj")
load("Robjs/tpm.mr.not_weighted.Robj")

tpm.mr <- tpm.mr[order(tpm.mr$rsqr, decreasing = T), ]
tpm.mr <- tpm.mr[order(tpm.mr$pval.best, decreasing = F), ]

head(data.frame(tpm.mr), n = 60)
# 
# sanity?
jgene <- "Ddc"

tpm.test <- subset(tpm.afe.avg, gene_name == jgene & tissue != "WFAT")
dat.mat <- dcast(tpm.test, tissue + amp + mean~ transcript_id, value.var = "tpm_norm.avg")
dat.proms <- subset(dat.mat, select = -c(tissue, amp, mean))
dat.H <- apply(dat.proms, 1, ShannonEntropy)

test <- GetPromoterUsage(tpm.test, jvar = "tpm_norm.avg", do.svd = F, append.tiss = TRUE, get.means = TRUE, get.entropy = TRUE)
test.svd <- GetPromoterUsage(tpm.test, do.svd = T, append.tiss = TRUE, get.means = TRUE)
CorrelateAmpPromMulti(tpm.test, jvar = "tpm.avg", do.svd = T, weighted = F, eps = 1e-10)
CorrelateAmpPromMulti(tpm.test, jvar = "tpm_norm.avg", do.svd = T, weighted = F, eps = 1e-10)

# cancor
DoCanCor(tpm.test, xvar = "tpm_norm.avg", yvar = c("amp"))
DoCanCor(tpm.test, xvar = "tpm_norm.avg", yvar = c("amp", "mean"))

m1 <- PlotGeneAcrossTissues(subset(dat.long, gene == jgene), convert.linear = F)
m2 <- ggplot(subset(tpm.afe.avg, gene_name == jgene), aes(x = transcript_id, y = tpm_norm.avg)) + geom_bar(stat = "identity") + facet_wrap(~tissue) + ggtitle(jgene)
m3 <- ggplot(subset(tpm.afe.avg, gene_name == jgene), aes(x = transcript_id, y = tpm.avg)) + geom_bar(stat = "identity") + facet_wrap(~tissue) + ggtitle(jgene)
jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
multiplot(m1, m2, m3, layout = jlayout)

plot(test$dat.mat.trans[, 2], test$dat.mat.trans[, 3])
text(test$dat.mat.trans[, 2], test$dat.mat.trans[, 3], labels = test$dat.mat.trans$tissue)

plot(test.svd$dat.mat.trans[, 2], test.svd$dat.mat.trans[, 3])
text(test.svd$dat.mat.trans[, 2], test.svd$dat.mat.trans[, 3], labels = test.svd$dat.mat.trans$tissue)

dat.mat <- dcast(tpm.test, tissue + amp + mean~ transcript_id, value.var = "tpm_norm.avg")
dat.mat.prom <- subset(dat.mat, select = -c(tissue, amp, mean))

dat.mat.prom <- sweep(dat.mat.prom, MARGIN = 1, STATS = rowMeans(dat.mat.prom), FUN = "-")
dat.mat.prom.s <- svd(dat.mat.prom)

dat.mat.prom.trans <- sweep(dat.mat.prom.s$u, MARGIN = 2, STATS = dat.mat.prom.s$d, FUN = "*")
plot(dat.mat.prom.s$u[, 1] * dat.mat.prom.s$d[1], dat.mat.prom.s$u[, 2]* dat.mat.prom.s$d[2])
text(dat.mat.prom.s$u[, 1] * dat.mat.prom.s$d[1], dat.mat.prom.s$u[, 2]* dat.mat.prom.s$d[2], labels = test.svd$dat.mat.trans$tissue)



# Cluster distance --------------------------------------------------------

for (jgene in head(tpm.mr, n = 20)$gene_name){
  tpm.test <- subset(tpm.afe.avg, gene_name == jgene & tissue != "WFAT")
  test <- GetPromoterUsage(tpm.test, jvar = "tpm_norm.avg", do.svd = F, append.tiss = TRUE, get.means = TRUE, get.entropy = TRUE)
  test.svd <- GetPromoterUsage(tpm.test, do.svd = T, append.tiss = TRUE, get.means = TRUE)
  
  dat <- subset(test.svd$dat.mat.trans, select = -tissue)
  dat.prom <- subset(dat, select = - amp)
  dat.amp <- dat$amp
  S <- 0
  
  for (i in seq(nrow(dat.prom))){
    prom1 <- dat.prom[i, ]
    amp1 <- dat.amp[i]
    for (j in i:nrow(dat.prom)){
      prom2 <- dat.prom[j, ]
      amp2 <- dat.amp[j]
      S <- S + PromAmpScore(prom1, prom2, amp1, amp2)
    }
  }
  print(jgene)
  print(S) 
}


# Fit Gaussians calculate likelihood --------------------------------------

jgene <- "Ddc"

RunFuzzyDistance(subset(tpm.afe.avg, gene_name == jgene & tissue != "WFAT"))

tpm.test <- subset(tpm.afe.avg, gene_name == jgene & tissue != "WFAT")
dat.mat <- dcast(tpm.test, tissue + amp + mean~ transcript_id, value.var = "tpm_norm.avg")
dat.proms <- subset(dat.mat, select = -c(tissue, amp, mean))
dat.H <- apply(dat.proms, 1, ShannonEntropy)

test <- GetPromoterUsage(tpm.test, jvar = "tpm_norm.avg", do.svd = F, append.tiss = TRUE, get.means = TRUE, get.entropy = TRUE)
test.svd <- GetPromoterUsage(tpm.test, do.svd = T, append.tiss = TRUE, get.means = TRUE)

plot(test$dat.mat.trans[, 2], test$dat.mat.trans[, 3])
text(test$dat.mat.trans[, 2], test$dat.mat.trans[, 3], labels = test$dat.mat.trans$tissue)

# proms <- test.svd$dat.mat.trans[, c(2, 3)]
proms <- subset(test.svd$dat.mat.trans, select = -c(amp, tissue))
amp <- test.svd$dat.mat.trans$amp
weights1 <- amp / max(amp)
weights2 <- 1 - weights1

# center1
mu1 <- colSums(sweep(proms, MARGIN = 1, STATS = weights1, FUN = "*")) / sum(weights1)

# center2
mu2 <- colSums(sweep(proms, MARGIN = 1, STATS = weights2, FUN = "*")) / sum(weights2)

plot(test.svd$dat.mat.trans[, 2], test.svd$dat.mat.trans[, 3])
text(test.svd$dat.mat.trans[, 2], test.svd$dat.mat.trans[, 3], labels = test.svd$dat.mat.trans$tissue)
points(c(mu1[1], mu2[1]), c(mu1[2], mu2[2]), pch = "*", col = "blue", cex = 5)

dist1 <- FuzzyDistance(proms, mu1, amp)
dist2 <- FuzzyDistance(proms, mu2, amp)
print(paste("Intracluster score", sum(dist1, dist2)))
print(paste("Interscore", sum((mu2 - mu1) ^ 2)))


# # Multivariate linear regression  -----------------------------------------
# 
# test <- subset(tpm.afe.avg, gene_name == jgene)
# test.mat <- dcast(test, tissue + amp + mean~ transcript_id, value.var = "tpm_norm.avg")
# test.mat.prom <- subset(test.mat, select = -c(tissue, amp, mean))
# test.mat.y <- subset(test.mat, select = c(tissue, amp, mean))
# 
# # reduce dim
# test.mat.prom <- sweep(test.mat.prom, MARGIN = 1, STATS = rowMeans(test.mat.prom), FUN = "-")
# test.mat.prom.s <- svd(test.mat.prom)
# eigvals <- test.mat.prom.s$d ^ 2 / sum(test.mat.prom.s$d ^ 2)
# eigvals.cum <- cumsum(eigvals)
# # threshold for number of eigvals to keep
# eigvals.sum.thres <- 0.9
# 
# 
# keep <- KeepUpToThres(eigvals.cum, eigvals.sum.thres, min.dim = 2)
# test.mat.prom.trans <- test.mat.prom.s$u[, keep] * test.mat.prom.s$d[keep]
# 
# plot(test.mat.prom.s$u[, 1], test.mat.prom.s$u[, 2])
# text(test.mat.prom.s$u[, 1], test.mat.prom.s$u[, 2], labels = test.mat$tissue)
# test.mat.trans <- data.frame(amp = test.mat$amp, test.mat.prom.trans)
# 
# fit.altprom <- lm(formula = amp ~ ., data = test.mat.trans, weights = (test.mat$mean))
# # fit.altprom <- lm(formula = amp ~ ., data = test.mat.trans)
# summary(fit.altprom)