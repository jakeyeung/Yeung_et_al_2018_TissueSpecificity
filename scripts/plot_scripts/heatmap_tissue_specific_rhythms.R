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

FindMaxAmp <- function(dat, amp.var = "mag.norm"){
  amp.max <- max(dat[[amp.var]])
  gene.max <- dat$gene[which(dat[[amp.var]] == amp.max)]
  return(data.frame(gene = gene.max, mag.norm = amp.max))
}


# Load --------------------------------------------------------------------

dat.long <- LoadArrayRnaSeq()

# filt.tiss <- c("Cere", "BS", "Hypo")
filt.tiss <- c()

# Fit relative amplitude --------------------------------------------------

# dat.fit <- FitRhythmicDatLong(dat.long)
load(file = "Robjs/dat.fit.Robj")

# Get amplitude relative to Bmal1 -----------------------------------------

ref.gene <- "Nr1d1"

# dat.fit.relamp <- GetRelamp(dat.fit, max.pval = 1e-3)
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)

dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))

# Get entropy measure -----------------------------------------------------


# 
# # entropy
# dat.entropy <- dat.fit.relamp %>%
#   #   subset(., min.pval <= 1e-5 & max.relamp > 0.1) %>%
#   subset(., max.exprs >= 4 & !tissue %in% filt.tiss) %>%
#   group_by(gene) %>%
#   #   summarise(entropy = ShannonEntropy(relamp / sum(relamp)))
#   summarise(entropy = ShannonEntropy(amp / sum(amp)))
# 
# plot(density(dat.entropy$entropy[which(!is.na(dat.entropy$entropy))]))
# 
# dat.entropy[order(dat.entropy$entropy), ]
# dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]
# 
# N <- floor(nrow(dat.entropy) / 3)
# n.genes <- nrow(dat.entropy)
# low.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy), ]$gene, n = N)
# high.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]$gene, n = N)
# med.entropy.genes <- dat.entropy[order(dat.entropy$entropy), ]$gene[N:(n.genes - N)]


# Get complex matrix ------------------------------------------------------

genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)

dat.complex <- TemporalToFrequencyDatLong(subset(dat.long, gene %in% genes.exprs), period = 24, n = 8, interval = 6, add.entropy.method = "array")

dat.complex$exprs.adj <- dat.complex$exprs.transformed * dat.complex$frac.weight

# adjust for "noise" and normalize to reference gene (Nr1d1)
ref.amps <- subset(dat.complex, gene == ref.gene, select = c(tissue, gene, exprs.adj))
ref.amps.dic <- hash(ref.amps$tissue, ref.amps$exprs.adj)

dat.complex$exprs.adj.norm <- mapply(function(tiss, exprs) exprs / ref.amps.dic[[tiss]], as.character(dat.complex$tissue), dat.complex$exprs.adj)
dat.complex$mod.exprs.adj <- Mod(dat.complex$exprs.adj)
dat.complex$mod.exprs.adj.norm <- Mod(dat.complex$exprs.adj.norm)

# Get entropy measure from adjusted magnitude -----------------------------



# ref gene should be max entropy as a check
dat.entropy <- dat.complex %>%
  subset(., gene %in% genes.exprs) %>%
  group_by(gene) %>%
  summarise(entropy = ShannonEntropy(Mod(exprs.adj.norm) / sum(Mod(exprs.adj.norm))))

plot(density(dat.entropy$entropy))

dat.entropy[order(dat.entropy$entropy), ]
dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]

N <- floor(nrow(dat.entropy) / 3)
n.genes <- nrow(dat.entropy)
low.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy), ]$gene, n = N)
high.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]$gene, n = N)
med.entropy.genes <- dat.entropy[order(dat.entropy$entropy), ]$gene[N:(n.genes - N)]

# Get entropy measure from Complex Matrix ---------------------------------

# # normalize by ref.gene
# dat.complex.sub <- subset(dat.complex, gene == ref.gene)
# ref.mag <- hash(as.character(dat.complex.sub$tissue), dat.sub$mag.adj)
# 
# dat.complex$exprs.norm.std <- mapply(function(tiss, exprs) exprs / ref.mag[[as.character(tiss)]], 
#                                     dat.complex$tissue, dat.complex$exprs.norm)
# 
# # normalize such that sum(Mod(x)) = 1
# dat.complex <- dat.complex %>%
#   group_by(gene) %>%
#   mutate(exprs.norm.std.norm = exprs.norm.std / sum(Mod(exprs.norm.std)))

# Get SVDs ----------------------------------------------------------------

s.all.norm <- SvdOnComplex(dat.complex, value.var = "exprs.adj")
s.low.norm <- SvdOnComplex(subset(dat.complex, gene %in% low.entropy.genes), value.var = "exprs.adj")
s.high.norm <- SvdOnComplex(subset(dat.complex, gene %in% high.entropy.genes), value.var = "exprs.adj")
s.med.norm <- SvdOnComplex(subset(dat.complex, gene %in% med.entropy.genes), value.var = "exprs.adj")

s.low <- SvdOnComplex(subset(dat.complex, gene %in% low.entropy.genes))
s.high <- SvdOnComplex(subset(dat.complex, gene %in% high.entropy.genes))
s.med <- SvdOnComplex(subset(dat.complex, gene %in% med.entropy.genes))
s.all <- SvdOnComplex(dat.complex)

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

# All
for (comp in seq(3)){
  eigens.all <- GetEigens(s.all, period = 24, comp = comp)
  multiplot(eigens.all$v.plot, eigens.all$u.plot, layout = jlayout)
}

# All.norm
for (comp in seq(3)){
  eigens.all.norm <- GetEigens(s.all.norm, period = 24, comp = comp)
  multiplot(eigens.all.norm$v.plot, eigens.all.norm$u.plot, layout = jlayout)
}

# High
for (comp in seq(3)){
  eigens.high <- GetEigens(s.high, period=24, comp=comp)
  multiplot(eigens.high$v.plot, eigens.high$u.plot, layout = jlayout)
}

# High.norm
for (comp in seq(2)){
  eigens.high.norm <- GetEigens(s.high.norm, period = 24, comp = comp)
  multiplot(eigens.high.norm$v.plot, eigens.high.norm$u.plot, layout = jlayout)
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

# Medium.norm
for (comp in seq(2)){
  eigens.med.norm <- GetEigens(s.med.norm, period=24, comp=comp)
  multiplot(eigens.med.norm$v.plot, eigens.med.norm$u.plot, layout = jlayout)
} 


# Can we discover new circadian genes? ------------------------------------

n.genes <- 200
# eg X6430573F11Rik
u1 <- sort(Mod(s.high.norm$u[, 1]), decreasing = TRUE)
u1.sub <- head(u1, n = n.genes)

dat.sub <- subset(dat.long, gene %in% names(u1.sub))

fpath <- file.path(outdir, "circadian_genes.pdf")
if (!file.exists(fpath)){
  pdf(fpath)
  barplot(u1.sub, names.arg = names(u1.sub), las = 2, horiz = FALSE, cex.names = 0.8)
  for (g in names(u1.sub)){
    print(PlotGeneAcrossTissues(subset(dat.sub, gene == g)))
  }
  dev.off()
}

sink(file = file.path(outdir, "circadian_genes.txt"))
for (g in names(u1.sub)){
  print(g)
}
sink()

# Can we identify tissue-specific rhythmic genes? -------------------------

n.genes <- 200
# eg X6430573F11Rik
for (comp in seq(3)){
  u1 <- sort(Mod(s.low.norm$u[, comp]), decreasing = TRUE)
  u1.sub <- head(u1, n = n.genes)
  
  dat.sub <- subset(dat.long, gene %in% names(u1.sub))
  outname <- paste0("tissue_spec_mod_", comp, "_", n.genes)
  fpath <- file.path(outdir, paste0(outname, ".pdf"))
  if (!file.exists(fpath)){
    pdf(fpath)
    barplot(u1.sub, names.arg = names(u1.sub), las = 2, horiz = FALSE, cex.names = 0.8)
    for (g in names(u1.sub)){
      print(PlotGeneAcrossTissues(subset(dat.sub, gene == g)))
    }
    dev.off()
    sink(file = file.path(outdir, paste0(outname, ".txt")))
    for (g in names(u1.sub)){
      print(g)
    }
    sink()
  }
}

# Heatmap low and high entropy genes -----------------------------------------------

pdf(file.path(outdir, "low_entropy_heatmap.nobrain.pdf"))
M <- LongToMat(subset(dat.complex, gene %in% low.entropy.genes & !tissue %in% filt.tiss), value.var = "mod.exprs.adj.norm")
# M <- LongToMat(subset(dat.fit.relamp, gene %in% low.entropy.genes & !tissue %in% filt.tiss))
PlotRelampHeatmap(M, paste0("Low entropy rhythmic genes. N=", length(low.entropy.genes)))
dev.off()

pdf(file.path(outdir, "high_entropy_heatmap.nobrain.pdf"))
# M <- LongToMat(subset(dat.fit.relamp, gene %in% high.entropy.genes & !tissue %in% filt.tiss))
M <- LongToMat(subset(dat.complex, gene %in% high.entropy.genes & !tissue %in% filt.tiss), value.var = "mod.exprs.adj.norm")
PlotRelampHeatmap(M, paste0("High entropy rhythmic genes. N=", length(high.entropy.genes)))
dev.off()

pdf(file.path(outdir, "medium_entropy_heatmap.nobrain.pdf"))
# M <- LongToMat(subset(dat.fit.relamp, gene %in% med.entropy.genes & !tissue %in% filt.tiss))
M <- LongToMat(subset(dat.complex, gene %in% med.entropy.genes & !tissue %in% filt.tiss), value.var = "mod.exprs.adj.norm")
PlotRelampHeatmap(M, paste0("Medium entropy rhythmic genes. N=", length(med.entropy.genes)))
dev.off()
# Plot top low entropy genes ---------------------------------------------

# heatmap
plot.n <- 500
outname <- paste0("low_entropy.", plot.n, ".pdf")
pdf(file.path(outdir, outname))
dat.sub <- subset(dat.long, gene %in% low.entropy.genes)
for (low.g in low.entropy.genes[1:plot.n]){
  print(low.g)
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == low.g)))
}
dev.off()


# Plot top high entropy genes ----------------------------------------------

outname <- paste0("high_entropy.", plot.n, ".pdf")
pdf(file.path(outdir, outname))
dat.sub <- subset(dat.long, gene %in% high.entropy.genes)
for (high.g in high.entropy.genes[1:plot.n]){
  print(high.g)
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == high.g)))
}
dev.off()


# Plot hierarchical clustering  -------------------------------------------

dat.complex$real <- Re(dat.complex$exprs.adj)
dat.complex$imag <- Im(dat.complex$exprs.adj)

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
