# 2016-07-01
# Jake Yeung
# Alternative promoter usage

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(hash)
library(mvtnorm)
library(reshape2)

source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/NcondsFunctions.R")

eps <- 1  # for log2 transform


# Function ----------------------------------------------------------------

comp <- function(x,mean1,sd1,mean2,sd2){
  # http://stats.stackexchange.com/questions/148385/overlap-between-two-normal-pdfs
  if (length(x)==1){ 
    outcome=min(dnorm(x,mean1,sd1),dnorm(x,mean2,sd2))}else{
      first=dnorm(x,mean1,sd1)
      second=dnorm(x,mean2,sd2)
      outcome=first*(first<second)+second*(first>=second)}
  return(outcome)}


PermutatePvalue <- function(tpm.test, n.perms, show.plot=FALSE){
  # http://stats.stackexchange.com/questions/30687/how-do-test-whether-two-multivariate-distributions-are-sampled-from-the-same-und
  set.seed(0)
  # TODO
  test.svd <- GetPromoterUsage(tpm.test, jvar = "tpm_norm.avg", do.svd = T, append.tiss = TRUE, get.means = TRUE, transcript_id = "transcript")
  proms <- subset(test.svd$dat.mat.trans, select = -c(amp, tissue))
  amp <- test.svd$dat.mat.trans$amp
  weights1 <- (amp - min(amp)) / (max(amp) - min(amp))
  weights2 <- 1 - weights1
  
  # calculate centers and dist
  mu1 <- colSums(sweep(proms, MARGIN = 1, STATS = weights1, FUN = "*")) / sum(weights1)
  mu2 <- colSums(sweep(proms, MARGIN = 1, STATS = weights2, FUN = "*")) / sum(weights2)
  jdist <- sum((mu2 - mu1) ^ 2)
  
  # PERMUTE
  jdists <- rep(NA, n.perms)
  for (i in seq(n.perms)){
    weights1.rand <- sample(weights1, size = length(weights1), replace = FALSE)
    weights2.rand <- 1 - weights1.rand
    mu1.rand <- colSums(sweep(proms, MARGIN = 1, STATS = weights1.rand, FUN = "*")) / sum(weights1.rand)
    mu2.rand <- colSums(sweep(proms, MARGIN = 1, STATS = weights2.rand, FUN = "*")) / sum(weights2.rand)
    jdist.rand <- sum((mu2.rand - mu1.rand) ^ 2)
    jdists[i] <- jdist.rand
  }
  
  if (show.plot){
    if (max(jdists) < jdist){
      plot(density(jdists), xlim=c(0, jdist * 1.2), main = jgene)
      abline(v = jdist)
    } else {
      plot(density(jdists), main = jgene)
      abline(v = jdist)
    }
  } 
  # calculate p-value
  N <- length(which(jdists > jdist))
  jpval <- N / length(jdists)
  return(jpval)
}






# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.bytranscript.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)

fits.orig <- fits.long.filt
fits.long.filt <- subset(fits.long.filt, method == "g=1001")

# Annotate fits
fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)
fits.long.filt$amp.avg <- mapply(GetAvgAmpFromParams, fits.long.filt$param.list, fits.long.filt$model)
fits.long.filt$phase.sd <- mapply(GetSdPhaseFromParams, fits.long.filt$param.list, fits.long.filt$model)
fits.long.filt$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.long.filt$param.list, fits.long.filt$model)
fits.long.filt$phase.avg <- mapply(GetPhaseFromParams, fits.long.filt$param.list, fits.long.filt$model)


dat.bytranscript <- StaggeredTimepointsLivKid(dat.bytranscript)
dat.bytranscript <- CollapseTissueGeno(dat.bytranscript)  # match fits.long.filt


dat.long <- StaggeredTimepointsLivKid(dat.long)

tpm.afe.avg <- dat.bytranscript %>%
  group_by(gene, tissue, time, geno) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1)) %>%
  group_by(gene, transcript, tissue, geno) %>%
  summarise(tpm_norm.avg = mean(tpm_normalized),
            tpm_norm.var = var(tpm_normalized),
            tpm.avg = mean(tpm),
            tpm.var = var(tpm))





# Positive control --------------------------------------------------------

jgene <- "Ddc"

jsub <- subset(dat.bytranscript, gene == jgene & geno == "SV129")

# PlotTpmAcrossTissues(jsub, jtitle = jgene, log2.transform = FALSE, transcript_id = "transcript")
# PlotGeneTissuesWTKO(subset(dat.long, gene == jgene))


# Count promoters ---------------------------------------------------------
tpm.counts <- subset(tpm.afe.avg, tissue == "Liver_SV129") %>%
  group_by(tissue, gene) %>%
  summarise(counts = length(transcript))
tpm.counts <- tpm.counts[order(tpm.counts$counts, decreasing = T), ]

counts.dic <- hash(as.character(tpm.counts$gene), tpm.counts$counts)
tpm.afe.avg$nprom <- sapply(as.character(tpm.afe.avg$gene), function(jgene) counts.dic[[jgene]])

# hash it
tiss <- c("Liver_SV129", "Kidney_SV129", "Liver_BmalKO", "Kidney_BmalKO")
keys.df <- fits.long.filt %>%
  group_by(gene) %>%
  do(GetGeneModelKeys(., tiss))
keys <- paste(keys.df$tissue, keys.df$gene, sep = ",")
vals <- as.numeric(keys.df$is.rhyth)
is.rhyth.dic <- hash(keys, vals)

# annotate it
tpm.afe.avg <- subset(tpm.afe.avg, gene %in% unique(fits.long.filt$gene) & geno == "SV129")

# 1 or 0 to signify whether it is rhythmic
tpm.afe.avg$amp <- mapply(function(jtiss, jgene) is.rhyth.dic[[paste(jtiss, jgene, sep = ",")]], as.character(tpm.afe.avg$tissue), as.character(tpm.afe.avg$gene))

# do only on subset of models??
jgenes <- unique(as.character(subset(fits.long.filt, model != "")$gene))
# jmods <- c("Liver_SV129", "Kidney_SV129")
# jgenes <- unique(as.character(subset(fits.long.filt, model %in% jmods)$gene))

# BEGIN: ANALYZE ON EACH TIME
tpm.afe.bysamp <- dat.bytranscript %>%
  group_by(gene, tissue, time, geno) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1)) %>%
  group_by(gene, transcript, tissue, geno, time) %>%
  summarise(tpm_norm.avg = mean(tpm_normalized),
            tpm_norm.var = var(tpm_normalized),
            tpm.avg = mean(tpm),
            tpm.var = var(tpm))
tpm.afe.bysamp <- subset(tpm.afe.bysamp, gene %in% unique(fits.long.filt$gene))
tpm.afe.bysamp$nprom <- sapply(as.character(tpm.afe.bysamp$gene), function(jgene) counts.dic[[jgene]])
tpm.afe.bysamp$amp <-  mapply(function(jtiss, jgene) is.rhyth.dic[[paste(jtiss, jgene, sep = ",")]], as.character(tpm.afe.bysamp$tissue), as.character(tpm.afe.bysamp$gene))

# Rename to just tissue, ignoring genotype
tpm.afe.bysamp$tissue.merge <- sapply(as.character(tpm.afe.bysamp$tissue), function(tiss) strsplit(tiss, "_")[[1]][[1]])

# amps of same tissue should be same (take maximum, if tissue in any genotype is rhythmic, then label it as rhythmic)
tpm.afe.bysamp <- tpm.afe.bysamp %>%
  group_by(gene, transcript, tissue.merge) %>%
  mutate(amp = max(amp))

# make tissue indivudal factors by adding geno and time
tpm.afe.bysamp$tissue <- paste(as.character(tpm.afe.bysamp$tissue), tpm.afe.bysamp$time, sep = "-")

tpm.gauss <- subset(tpm.afe.bysamp, nprom > 1 & gene %in% jgenes) %>%
  group_by(gene) %>%
  do(sigs = CalculateGaussianCenters(., transcript_id = "transcript"))

tpm.gauss2 <- subset(tpm.gauss, !is.na(sigs)) %>%
  group_by(gene) %>%
  do(CalculateGaussianDists(.))

# END: ANALYZE ON EACH TIME

genelist <- tpm.gauss2$gene
tpm.gauss <- cbind(as.data.frame(tpm.gauss2), as.data.frame(subset(tpm.gauss, gene %in% genelist, select = -gene)))

tpm.gauss <- tpm.gauss[order(tpm.gauss$center.dists, decreasing = TRUE), ]

head(data.frame(subset(tpm.gauss, select=-sigs)), n = 50)



# Testing single gene, commented out
# jgene <- "Mreg"
# # tpm.afe.bysamp <- subset(tpm.afe.bysamp, gene %in% unique(fits.long.filt$gene) & geno == "SV129")
# jtest <- subset(tpm.afe.bysamp, gene == jgene)
# jtest$tissue <- paste(jtest$tissue, jtest$time, sep = "-")
# out <- GetPromoterUsage(jtest, transcript_id = "transcript")
# 
# proms <- out$dat.mat.trans
# jtitle <- jgene
# plot(proms[, 2], proms[, 3], main = jtitle, xlab = "Promoter usage (1st component)", ylab = "Promoter usage (2nd component)")
# par(cex.axis=1, cex.lab=2, cex.main=2, cex.sub=1)
# 
# jcols <- vector(length = nrow(proms))
# jcols[which(proms$amp == 1)] <- "blue"
# jcols[which(proms$amp == 0)] <- "red"
# text(proms[, 2], proms[, 3], labels = proms$tissue, col = jcols)
# draw ellipse
# if (draw.ellipse){
#   try(lines(ellipse(tpm.gauss.sigs$sig1$cov, level = 0.5, centre = tpm.gauss.sigs$sig1$center), type='l', col = "blue"), silent = T)
#   try(lines(ellipse(tpm.gauss.sigs$sig2$cov, level = 0.5, centre = tpm.gauss.sigs$sig2$center), type='l', col = "red"), silent = T)
# }

# SANITY ------------------------------------------------------------------

# filter for similarly expressed genes (both high)

dat.min <- subset(dat.long, geno == "SV129") %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs)) %>%
  group_by(gene) %>%
  summarise(exprs.min = min(exprs.mean))

jcutoff <- 2
jgenes.min <- as.character(subset(dat.min, exprs.min >= jcutoff)$gene)

head(data.frame(subset(tpm.gauss, gene %in% jgenes.min, select=-sigs)), n = 50)
# 
# # TESTING
# jgene <- "Slc45a3"
# jgene <- "Gm14327"
# jgene <- "Pnrc2"
# PromoterSpacePlots(tpm.afe.bysamp, jgene = jgene, transcript_id = "transcript", use.weights = FALSE)
# PlotGeneTissuesWTKO(subset(dat.long, gene == jgene), jtitle = jgene)
# 
# jgene <- "Npnt"
# jgene <- "Osgin1"
# jgene <- "Zranb1"
# jgene <- "Aox3"
# jgene <- "Upp2"
# jgene <- "9530068E07Rik"
# jgene <- "Prkd3"
# jgene <- "Sh3bgrl2"
# jgene <- "Por"
# jgene <- "Aqp9"
# jgene <- "Npnt"
# jgene <- "Wdtc1"
# jgene <- "Zranb1"
# jgene <- "Psen2"
# jgene <- "Prkd3"
# jgene <- "Insig2"
# jgene <- "Upp2"
# jgene <- "Ddc"
# jgene <- "Usp2"
# jgene <- "Slc45a3"

genes <- as.character(head(data.frame(subset(tpm.gauss, gene %in% jgenes, select=-sigs)), n = 50)$gene)


# pdf("plots/alternative_exon_usage/liver_kidney.atger_nestle.bugfixed.bytime.pdf")
# for (jgene in genes){
#   print(jgene)
#   jmodel <- as.character(subset(fits.long.filt, gene == jgene)$model[[1]])
#   tx <- GetPromoterUsage(subset(tpm.afe.bysamp, gene == jgene), transcript_id = "transcript", get.entropy=FALSE, return.transcripts=TRUE)[[1]]$transcript
#   print(PlotGeneTissuesWTKO(subset(dat.long, gene == jgene)) + ggtitle(paste(jgene, jmodel)))
#   PromoterSpacePlots(tpm.afe.bysamp, jgene = jgene, transcript_id = "transcript", use.weights = FALSE)
#   print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & transcript %in% tx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
#   print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
# }
# dev.off()

# PromoterSpacePlots.nostics(subset(tpm.gauss, gene == jgene)$sigs[[1]], jgene, draw.ellipse = F)
# subset(fits.orig, gene == jgene)


# Calculate p-value for each gene: genome wide ----------------------------

genes.filt <- unique(as.character(subset(tpm.gauss, gene %in% jgenes)$gene))

tpm.afe.bysamp <- subset(tpm.afe.bysamp, gene %in% genes.filt) %>%
  group_by(gene) %>%
  do(pval = PermutatePvalue(., 10000, show.plot=FALSE))

save(tpm.gauss, tpm.afe.bysamp, file = "Robjs/liver_kidney_atger_nestle/alt_promoter_analysis_by_time.Robj")

# Calculate p-value from each gene ----------------------------------------

# jgene <- "Gm14327"
# jgene <- "Spns2"  # 50th
# jgene <- "Ube2g2"  # 5000th
# jgene <- "Fbxw9"  # 1000th
# jgene <- "Nasp"  # 500th
# jgene <- "Srpk1"  # 250th
# 
# PromoterSpacePlots(tpm.afe.bysamp, jgene = jgene, transcript_id = "transcript", use.weights = FALSE)
# tpm.test <- subset(tpm.afe.bysamp, gene == "Slc45a3")
# PermutatePvalue(tpm.test, 10000, show.plot=FALSE)

