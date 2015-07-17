# 2015-07-15
# Jake Yeung
# Show heatmap of tissue-specific rhythms without being dependent on p-values

setwd("~/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/entropy_genes"

start.time <- Sys.time()

library(dplyr)
library(gplots)
library(reshape2)

# Functions ---------------------------------------------------------------
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LongToMat.R")
source("scripts/functions/ShannonEntropy.R")
source("scripts/functions/DataHandlingFunctions.R")
clockgenes <- GetClockGenes()

# Load --------------------------------------------------------------------

dat.long <- LoadArrayRnaSeq()

# Fit relative amplitude --------------------------------------------------

dat.fit <- FitRhythmicDatLong(dat.long)

# Get amplitude relative to Bmal1 -----------------------------------------

ref.gene <- "Arntl"

dat.fit.relamp <- GetRelamp(dat.fit, max.pval = 1e-3)
# dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)

dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))

# Get entropy measure -----------------------------------------------------

# filt.tiss <- c("Cere", "BS", "Hypo")
filt.tiss <- c()

# entropy
dat.entropy <- dat.fit.relamp %>%
#   subset(., min.pval <= 1e-5 & max.relamp > 0.1) %>%
  subset(., max.exprs >= 4 & max.relamp >= 0.1 & min.pval <= 1e-3 & !tissue %in% filt.tiss) %>%
  group_by(gene) %>%
  summarise(entropy = ShannonEntropy(relamp / sum(relamp)))

plot(density(dat.entropy$entropy[which(!is.na(dat.entropy$entropy))]))

dat.entropy[order(dat.entropy$entropy), ]
dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]

N <- 2000
n.genes <- nrow(dat.entropy)
low.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy), ]$gene, n = N)
high.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]$gene, n = N)
med.entropy.genes <- dat.entropy[order(dat.entropy$entropy), ]$gene[N:(n.genes - N)]


# Get complex matrix ------------------------------------------------------

# ~ a minute
# dat.complex <- TemporalToFrequency(subset(dat.long, gene %in% dat.entropy$gene))

# take lowest sampling rate n=8, interval=6 -> rnaseq
# slow ~ 5 minutes
# dat.complex <- dat.long %>%
#   subset(., gene %in% dat.entropy$gene) %>%
#   group_by(gene, tissue) %>%
#   do(TemporalToFrequency2(., period = 24, n = 8, interval = 6, add.entropy = TRUE))

# ~ 1 minute
# dat.complex <- dat.long %>%
#   subset(., gene %in% dat.entropy$gene) %>%
#   group_by(gene, tissue) %>%
#   do(TemporalToFrequency2(.))

dat.complex <- TemporalToFrequencyDatLong(subset(dat.long, gene %in% dat.entropy$gene), period = 24, n = 8, interval = 6, add.entropy.method = "array")
dat.complex$magnitude <- Mod(dat.complex$exprs.transformed)
dat.complex$exprs.norm <- dat.complex$exprs.transformed * dat.complex$log.H.weight

s.all.norm <- SvdOnComplex(dat.complex, value.var = "exprs.norm")
s.low.norm <- SvdOnComplex(subset(dat.complex, gene %in% low.entropy.genes), value.var = "exprs.norm")

s.low <- SvdOnComplex(subset(dat.complex, gene %in% low.entropy.genes))
s.high <- SvdOnComplex(subset(dat.complex, gene %in% high.entropy.genes))
s.med <- SvdOnComplex(subset(dat.complex, gene %in% med.entropy.genes))
s.all <- SvdOnComplex(dat.complex)

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

# All
for (comp in seq(2)){
  eigens.all <- GetEigens(s.all, period = 24, comp = comp)
  multiplot(eigens.all$v.plot, eigens.all$u.plot, layout = jlayout)
}

# All.norm
for (comp in seq(2)){
  eigens.all.norm <- GetEigens(s.all.norm, period = 24, comp = comp)
  multiplot(eigens.all.norm$v.plot, eigens.all.norm$u.plot, layout = jlayout)
}

# High
for (comp in seq(1)){
  eigens.high <- GetEigens(s.high, period=24, comp=comp)
  multiplot(eigens.high$v.plot, eigens.high$u.plot, layout = jlayout)
}

# Low
for (comp in seq(3)){
  eigens.low <- GetEigens(s.low, period=24, comp=comp)
  multiplot(eigens.low$v.plot, eigens.low$u.plot, layout = jlayout)
}

# Low.norm
for (comp in seq(3)){
  eigens.low.norm <- GetEigens(s.low.norm, period=24, comp=comp)
  multiplot(eigens.low.norm$v.plot, eigens.low.norm$u.plot, layout = jlayout)
}  

# Medium
for (comp in seq(2)){
  eigens.med <- GetEigens(s.med, period = 24, comp = comp)
  multiplot(eigens.med$v.plot, eigens.med$u.plot, layout = jlayout)
}

# Heatmap low and high entropy genes -----------------------------------------------

pdf(file.path(outdir, "low_entropy_heatmap.nobrain.pdf"))
M <- LongToMat(subset(dat.fit.relamp, gene %in% low.entropy.genes & !tissue %in% filt.tiss))
PlotRelampHeatmap(M, paste0("Low entropy rhythmic genes. N=", length(low.entropy.genes)))
dev.off()

pdf(file.path(outdir, "high_entropy_heatmap.nobrain.pdf"))
M <- LongToMat(subset(dat.fit.relamp, gene %in% high.entropy.genes & !tissue %in% filt.tiss))
PlotRelampHeatmap(M, paste0("High entropy rhythmic genes. N=", length(high.entropy.genes)))
dev.off()

pdf(file.path(outdir, "medium_entropy_heatmap.nobrain.pdf"))
M <- LongToMat(subset(dat.fit.relamp, gene %in% med.entropy.genes & !tissue %in% filt.tiss))
PlotRelampHeatmap(M, paste0("Medium entropy rhythmic genes. N=", length(med.entropy.genes)))
dev.off()
# Plot top low entropy genes ---------------------------------------------

# heatmap
plot.n <- 500
outname <- paste0("low_entropy.", plot.n, ".nobrain.pdf")
pdf(file.path(outdir, outname))
dat.sub <- subset(dat.long, gene %in% low.entropy.genes)
for (low.g in low.entropy.genes[1:plot.n]){
  print(low.g)
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == low.g)))
}
dev.off()


# Plot top high entropy genes ----------------------------------------------

outname <- paste0("high_entropy.", plot.n, ".nobrain.pdf")
pdf(file.path(outdir, outname))
dat.sub <- subset(dat.long, gene %in% high.entropy.genes)
for (high.g in high.entropy.genes[1:plot.n]){
  print(high.g)
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == high.g)))
}
dev.off()


# Plot hierarchical clustering  -------------------------------------------

dat.complex$real <- Re(dat.complex$exprs.transformed)
dat.complex$imag <- Im(dat.complex$exprs.transformed)

M.real <- t(LongToMat(dat.complex, value.var = "real"))
M.imag <- t(LongToMat(dat.complex, value.var = "imag"))
M.all <- cbind(M.real, M.imag)
hc <- hclust(dist(M.all, method = "euclidean"))
plot(hc, main = "Clustering on fourier component T=24h")

# Cluster on normalized ---------------------------------------------------

# dat.ref <- subset(dat.complex, gene == ref.gene)
# dic <- hash(as.character(dat.ref$tissue), dat.ref$magnitude)

dat.ref <- subset(dat.fit.relamp, relamp == 1)
dic <- hash(as.character(dat.ref$tissue), dat.ref$amp / 2)

dat.complex.norm <- dat.complex %>%
  group_by(tissue) %>%
  do(NormalizeComplexMat(., dic))

dat.complex.norm$real <- Re(dat.complex.norm$exprs.norm)
dat.complex.norm$imag <- Im(dat.complex.norm$exprs.norm)

M.real <- t(LongToMat(dat.complex.norm, value.var = "real"))
M.imag <- t(LongToMat(dat.complex.norm, value.var = "imag"))
M.all <- cbind(M.real, M.imag)
hc <- hclust(dist(M.all, method = "euclidean"))
plot(hc, main = "Clustering on fourier component T=24h")

# Done --------------------------------------------------------------------

print(Sys.time() - start.time)
