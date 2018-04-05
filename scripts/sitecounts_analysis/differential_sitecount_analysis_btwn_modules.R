# Differential sitecount analysis: different within tissue-specific modules
# 2015-11-20

setwd("~/projects/tissue-specificity/")

library(hash)
library(dplyr)
library(ggplot2)
library(reshape2)

source("scripts/functions/LoadSitecounts.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/DifferentialSitecountsFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")

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


load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj", verbose=T)


# Normalize sum -----------------------------------------------------------

N.sum <- N %>%
  group_by(tissue) %>%
  summarise(sitecount.sum = sum(motevo.value))
ggplot(N.sum, aes(y = sitecount.sum, x = tissue)) + geom_boxplot()  #liver has in general higher sitecounts
sum.dic <- hash(as.character(N.sum$tissue), N.sum$sitecount.sum)

# normalize across factors and genes
N$sitecount.norm <- mapply(function(jtiss, jval){
  return(10 ^ 6 * jval / sum.dic[[jtiss]])
}, as.character(N$tissue), N$motevo.value)

# normalize across genes
N.sum.gene <- N %>%
  group_by(tissue, motif) %>%
  summarise(sitecount.sum = sum(sitecount.norm))
key <- paste(as.character(N.sum.gene$tissue), as.character(N.sum.gene$motif), sep = ";")
val <- N.sum.gene$sitecount.sum
sum.gene.dic <- hash(key, val)

N$sitecount.norm.bymotif <- mapply(function(jtiss, jmotif, jval){
  return(10 ^ 3 * jval / sum.gene.dic[[paste(jtiss, jmotif, sep = ";")]])
}, as.character(N$tissue), as.character(N$motif), N$sitecount.norm)

# Do 2 by 2 tables ---------------------------------------------------------

genes.liver <- as.character(subset(fits.best, model == "Liver")$gene)
genes.bfat <- as.character(subset(fits.best, model == "BFAT")$gene)
genes.flat <- as.character(subset(fits.best, model == "")$gene)

jtiss <- "Liver"
jmodel <- "Liver"
genes.rhyth <- as.character(subset(fits.best, model == jmodel)$gene)

N.sub <- subset(N, gene %in% genes.rhyth)

motifs <- as.character(unique(N.sub$motif))

N.sub$sitecount <- N.sub$motevo.value
N.sub$model <- sapply(as.character(N.sub$tissue), function(x){
  if (x == jtiss){
    return(jtiss)
  } else {
    return("Flat")
  }
})

test.out <- N.sub %>%
  group_by(motif) %>%
  do(FisherTestSitecounts(., cutoff = 1.1))
test.out[order(test.out$p.value), ]


start <- Sys.time()
cutoffs <- seq(from = 0.05, to = 0.2, by = 0.05)
test.all <- data.frame()
for (cutoff in cutoffs){
  print(cutoff)
  test.out <- N.sub %>%
    group_by(motif) %>%
    do(FisherTestSitecounts(., cutoff, sitecount.col = "sitecount.norm", show.table = FALSE))
  test.out$cutoff <- cutoff
  test.all <- rbind(test.all, test.out)
}
print(Sys.time() - start)

test.sum <- test.all %>%
  group_by(motif) %>%
  summarise(odds.ratio = mean(odds.ratio), p.value.minuslog = mean(-log10(p.value)))

ggplot(test.sum, aes(y = p.value.minuslog, x = odds.ratio, label = motif)) + geom_point() + geom_text()

FisherTestSitecounts(subset(N.sub, motif == "CUX2.p2"), cutoff=1, show.table = TRUE)
FisherTestSitecounts(subset(N.sub, motif == "bHLH_family.p2"), cutoff=1, show.table = TRUE)
FisherTestSitecounts(subset(N.sub, motif == "HNF1A.p2"), cutoff=1, show.table = TRUE)
FisherTestSitecounts(subset(N.sub, motif == "RXRG_dimer.p3"), cutoff=1, show.table = TRUE)
FisherTestSitecounts(subset(N.sub, motif == "AIRE.p2"), cutoff=1, show.table = TRUE)


# Reduce clutter ----------------------------------------------------------

oddsrat.cutoff <- 1.7
test.sum$label <- mapply(function(jmotif, oddsrat){
  if (oddsrat > oddsrat.cutoff){
    return(as.character(jmotif))
  } else {
    return("")
  }
}, test.sum$motif, test.sum$odds.ratio)

ggplot(test.sum, aes(y = p.value.minuslog, x = odds.ratio, label = label)) + geom_point() + geom_text(position = position_jitter(w = 0.1, h = 0.1))

library(wordcloud)
textplot(x = test.sum$odds.ratio, y = test.sum$p.value.minuslog, words = test.sum$label, xlab = "Odds ratio", ylab = "-log10(pvalue)")
points(x = test.sum$odds.ratio, y = test.sum$p.value.minuslog, pch = ".", cex = 3)


# Enrichment for pairs of motifs? -----------------------------------------

# RXRG and RORA (and bHLH?) cross with all others
sitecount.col <- "sitecount.norm"
N.rora <- subset(N.sub, motif == "RORA.p2")
key.rora <- paste(as.character(N.rora$gene), as.character(N.rora$tissue), sep = ";")
val.rora <- N.rora[[sitecount.col]]
dic.rora <- hash(key.rora, val.rora)

N.rxrg <- subset(N.sub, motif == "RXRG_dimer.p3")
key.rxrg <- paste(as.character(N.rxrg$gene), as.character(N.rxrg$tissue), sep = ";")
val.rxrg <- N.rxrg[[sitecount.col]]
dic.rxrg <- hash(key.rxrg, val.rxrg)

N.cux2 <- subset(N.sub, motif == "CUX2.p2")
key.cux2 <- paste(as.character(N.cux2$gene), as.character(N.cux2$tissue), sep = ";")
val.cux2 <- N.cux2[[sitecount.col]]
dic.cux2 <- hash(key.cux2, val.cux2)

N.sub$sitecount.rora <- mapply(function(jgene, jtiss){
  return(dic.rora[[paste(as.character(jgene), as.character(jtiss), sep = ";")]])
}, N.sub$gene, N.sub$tissue)

N.sub$sitecount.rxrg <- mapply(function(jgene, jtiss){
  return(dic.rxrg[[paste(as.character(jgene), as.character(jtiss), sep = ";")]])
}, N.sub$gene, N.sub$tissue)

N.sub$sitecount.cux2 <- mapply(function(jgene, jtiss){
  return(dic.cux2[[paste(as.character(jgene), as.character(jtiss), sep = ";")]])
}, N.sub$gene, N.sub$tissue)

N.sub$sitecount.roraXmotif <- N.sub$sitecount * N.sub$sitecount.rora
N.sub$sitecount.rxrgXmotif <- N.sub$sitecount * N.sub$sitecount.rxrg
N.sub$sitecount.cux2Xmotif <- N.sub$sitecount * N.sub$sitecount.cux2

jmotif <- "HNF4A_NR2F1.2.p2"
FisherTestSitecounts(subset(N.sub, motif == jmotif), cutoff=1.5, sitecount.col = "sitecount.roraXmotif", show.table = TRUE)
FisherTestSitecounts(subset(N.sub, motif == jmotif), cutoff=1.5, sitecount.col = "sitecount.rxrgXmotif", show.table = TRUE)
FisherTestSitecounts(subset(N.sub, motif == jmotif), cutoff=1.5, sitecount.col = "sitecount.cux2Xmotif", show.table = TRUE)
FisherTestSitecounts(subset(N.sub, motif == jmotif), cutoff=1.5, sitecount.col = "sitecount", show.table = TRUE)


start <- Sys.time()
sitecount.col <- "sitecount.rxrgXmotif"
cutoffs <- seq(from = 0.5, to = 2, by = 0.5)
test.all.cross <- data.frame()
for (cutoff in cutoffs){
  print(cutoff)
  test.out.cross <- N.sub %>%
    group_by(motif) %>%
    do(FisherTestSitecounts(., cutoff, sitecount.col = sitecount.col, show.table = FALSE))
  test.out.cross$cutoff <- cutoff
  test.all.cross <- rbind(test.all.cross, test.out.cross)
}
print(Sys.time() - start)

test.sum.cross <- test.all.cross %>%
  group_by(motif) %>%
  summarise(odds.ratio = mean(odds.ratio), p.value.minuslog = mean(-log10(p.value)))

ggplot(test.sum.cross, aes(y = p.value.minuslog, x = odds.ratio, label = motif)) + geom_point() + geom_text()

jmotif <- "CEBPA.B_DDIT3.p2"
jmotif <- "RXRG_dimer.p3"
FisherTestSitecounts(subset(N.sub, motif == jmotif), cutoff=1.5, sitecount.col = "sitecount", show.table = TRUE)
FisherTestSitecounts(subset(N.sub, motif == jmotif), cutoff=1.5, sitecount.col = "sitecount.cux2Xmotif", show.table = TRUE)

