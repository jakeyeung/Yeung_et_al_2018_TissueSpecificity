# 2015-07-15
# Jake Yeung
# Show heatmap of tissue-specific rhythms without being dependent on p-values
# instead of entropy use variance

setwd("~/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/entropy_genes"

# HARDCODED USER-MODIFIABLE\
writedir.prefix <- "/home/yeung/projects/tissue-specificity/data/gene_lists/entropy_genes_split_by_"
filt.tiss <- c("Liver", "BFAT", "Mus")
# filt.tiss <- c("Liver")
# filt.tiss <- c("Liver", "Adr", "Aorta", "BFAT", "BS", "Hypo", "WFAT")
# filt.tiss <- c()
N.groups <- 2
if (length(filt.tiss) == 0){
  tissues.filtered <- "all_tissues"
} else {
  tissues.filtered <- paste0("no", paste(filt.tiss, collapse=""))
}
writedir <- paste0(writedir.prefix, N.groups, ".", tissues.filtered)
print(paste("Writing to:", writedir))
# writedir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/entropy_genes_split_by_3.noliver"

dir.create(writedir)

start.time <- Sys.time()

library(dplyr)
library(gplots)
library(reshape2)

# Functions ---------------------------------------------------------------
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
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

GetHighLowEntropy <- function(dat.byvar, split.by = 2, colname = "entropy"){
  N <- floor(nrow(dat.byvar) / split.by)
  n.genes <- nrow(dat.byvar)
  low.entropy.genes <- head(dat.byvar[order(dat.byvar[[colname]]), ]$gene, n = N)
  high.entropy.genes <- head(dat.byvar[order(dat.byvar[[colname]], decreasing = TRUE), ]$gene, n = N)
  med.entropy.genes <- dat.byvar[order(dat.byvar[[colname]]), ]$gene[N:(n.genes - N)]
  return(list(lowE=low.entropy.genes, medE=med.entropy.genes, highE=high.entropy.genes))
}

# Load --------------------------------------------------------------------

dat.long <- LoadArrayRnaSeq(fix.rik.xgene = TRUE)

# filt.tiss <- c("Cere", "BS", "Hypo")

# Fit relative amplitude --------------------------------------------------

# dat.fit <- FitRhythmicDatLong(dat.long)
load(file = "Robjs/dat.fit.Robj")

dat.fit <- subset(dat.fit, ! tissue %in% filt.tiss)

# Get amplitude relative to Bmal1 -----------------------------------------

ref.gene <- "Nr1d1"

# dat.fit.relamp <- GetRelamp(dat.fit, max.pval = 1e-3)
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)

dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))


# Get complex matrix ------------------------------------------------------

genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)

# dat.complex <- TemporalToFrequencyDatLong(subset(dat.long, gene %in% genes.exprs), period = 24, n = 8, interval = 6, add.entropy.method = "array")
# save(dat.complex, file = "Robjs/dat.complex.raw.Robj")
load(file = "Robjs/dat.complex.raw.Robj", verbose = TRUE)

dat.complex <- subset(dat.complex, ! tissue %in% filt.tiss)

dat.complex$exprs.adj <- dat.complex$exprs.transformed * dat.complex$frac.weight
dat.complex$mod.exprs.adj <- Mod(dat.complex$exprs.adj)

# adjust for "noise" and normalize to reference gene (Nr1d1)
ref.amps <- subset(dat.complex, gene == ref.gene, select = c(tissue, gene, exprs.adj))
ref.amps.dic <- hash(ref.amps$tissue, Mod(ref.amps$exprs.adj))

dat.complex$exprs.adj.norm <- mapply(function(tiss, exprs) exprs / ref.amps.dic[[tiss]], as.character(dat.complex$tissue), dat.complex$exprs.adj)
dat.complex$mod.exprs.adj.norm <- Mod(dat.complex$exprs.adj.norm)

# # adjust to total noise. Commented out because I like Nr1d1 to have mod.exprs.adj.norm = 1
# ref.amps.total <- dat.complex %>%
#   group_by(tissue) %>%
#   summarise(total_amp = sum(Mod(exprs.adj) ^ 2))
# ref.amps.total.dic <- hash(ref.amps.total$tissue, ref.amps.total$total_amp)
# dat.complex$exprs.adj.norm <- mapply(function(tiss, exprs) exprs / ref.amps.total.dic[[tiss]], as.character(dat.complex$tissue), dat.complex$exprs.adj)
# dat.complex$mod.exprs.adj.norm <- Mod(dat.complex$exprs.adj.norm)

# Find tissue-specific and tissue-wide genes by variance ------------------

dat.var <- dat.complex %>%
  group_by(gene) %>%
  summarise(var.noise.norm = var(mod.exprs.adj / sum(mod.exprs.adj)),
            var.noise.norm.log10 = -log10(var.noise.norm))

# dat.var <- dat.complex %>%
#   group_by(gene) %>%
#   summarise(var.noise.norm = var(mod.exprs.adj),
#             var.noise.norm.log10 = -log10(var.noise.norm))

plot(density(dat.var$var.noise.norm.log10),
     main = "Distribution of -log10(var)", cex.main = 1.75, cex.lab = 1.75, cex.axis = 1.75,
     xlab = "-log10(var)")
abline(v = sort(dat.var$var.noise.norm.log10)[5244])
abline(v = sort(dat.var$var.noise.norm.log10)[length(dat.var$var.noise.norm.log10) - 5244])

dat.var[order(dat.var$var.noise.norm.log10, decreasing = TRUE), ]


# Get entropy measure from adjusted magnitude -----------------------------
# normalized by Nr1d1

# ref gene should be max entropy as a check
dat.entropy <- dat.complex %>%
  subset(., gene %in% genes.exprs) %>%
  group_by(gene) %>%
  summarise(entropy = ShannonEntropy(Mod(exprs.adj.norm) / sum(Mod(exprs.adj.norm))))

plot(density(dat.entropy$entropy))

dat.entropy[order(dat.entropy$entropy), ]
dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]

# By entropy: not normalized by Nr1d1 -------------------------------------

dat.entropy.notadj <- dat.complex %>%
  subset(., gene %in% genes.exprs) %>%
  group_by(gene) %>%
  summarise(entropy = ShannonEntropy(Mod(exprs.adj) / sum(Mod(exprs.adj))))

plot(density(dat.entropy.notadj$entropy))

dat.entropy.notadj[order(dat.entropy.notadj$entropy), ]
dat.entropy.notadj[order(dat.entropy.notadj$entropy, decreasing = TRUE), ]

# Get low, med, high entropy genes ----------------------------------------

entropy.genes <- GetHighLowEntropy(dat.entropy, split.by = N.groups, colname = "entropy")
entropy.genes.notadj <- GetHighLowEntropy(dat.entropy.notadj, split.by = N.groups, colname = "entropy")
entropy.genes.var <- GetHighLowEntropy(dat.var, split.by = N.groups, colname = "var.noise.norm.log10")

low.entropy.genes.all <- list(entropy = entropy.genes$lowE, entropynotadj = entropy.genes.notadj$lowE, Var = entropy.genes.var$lowE)
high.entropy.genes.all <- list(entropy = entropy.genes$highE, entropynotadj = entropy.genes.notadj$highE, Var = entropy.genes.var$highE)

# Get SVDs ----------------------------------------------------------------

s.low.norm.all <- lapply(low.entropy.genes.all, function(g) SvdOnComplex(subset(dat.complex, gene %in% g), value.var = "exprs.adj"))
s.high.norm.all <- lapply(high.entropy.genes.all, function(g) SvdOnComplex(subset(dat.complex, gene %in% g), value.var = "exprs.adj"))


# SVD on low and high entropy genes ---------------------------------------

# HIGH

lapply(s.high.norm.all, PlotFirstNComponents, comps = 1)
lapply(s.low.norm.all, PlotFirstNComponents, comps = 1)


# Write gene lists to file to run MARA ------------------------------------


# for (E in names(low.entropy.genes.all)){
#   outname <- paste0(E, ".low_entropy.gene_list")
#   sink(file = file.path(writedir, outname))
#   for (g in low.entropy.genes.all[[E]]){
#     cat(g)
#     cat("\n")
#   }
#   sink()
#   outname <- paste0(E, ".high_entropy.gene_list")  # assume E is common between low and high
#   sink(file = file.path(writedir, outname))
#   for (g in high.entropy.genes.all[[E]]){
#     cat(g)
#     cat("\n")
#   }
#   sink()
# }

# Var only
E = "Var"
outname <- paste0(E, ".low_entropy.gene_list")
sink(file = file.path(writedir, outname))
for (g in low.entropy.genes.all[[E]]){
  cat(g)
  cat("\n")
}
sink()
outname <- paste0(E, ".high_entropy.gene_list")  # assume E is common between low and high
sink(file = file.path(writedir, outname))
for (g in high.entropy.genes.all[[E]]){
  cat(g)
  cat("\n")
}
sink()
outname <- paste0("expressed_genes.gene_list")  # assume E is common between low and high
sink(file = file.path(writedir, outname))
for (g in genes.exprs){
  cat(g)
  cat("\n")
}
sink()

# Run MARA on these sets of genes -----------------------------------------


