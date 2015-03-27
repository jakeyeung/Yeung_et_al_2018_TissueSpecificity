# Jake Yeung
# January 7 2015
# To eliminate genes with "noise", which shows up in projection and subsequent svd analysis, let's do projection only
# on genes that show rhythmicity in at least one tissue.
# 
# Dependencies: output of find.oscillating.genes.R is used to detect rhythmic genes.

library(plyr)
library(reshape2)
library(doMC)
library(PhaseHSV)
library(ggplot2)


# Parameters to play ------------------------------------------------------

## Specify the number of cores
mc <- 30

# Minimum log2 expression to be considered "expressed" in RNA-Seq
cutoff <- 5


clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')

clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')

# My Functions ------------------------------------------------------------

scripts.dir <- "scripts"
funcs.dir <- "functions"
source(file.path(scripts.dir, funcs.dir, "MergeAndProject.R"))
source(file.path(scripts.dir, funcs.dir, "GetTissueTimes.R"))
source(file.path(scripts.dir, funcs.dir, "FourierFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "LoadAndHandleData.R"))
source(file.path(scripts.dir, funcs.dir, "DataHandlingFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "PlotFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "GetTopNValues.R"))
source(file.path(scripts.dir, funcs.dir, "OuterComplex.R"))
source(file.path(scripts.dir, funcs.dir, "OrderPhaseMatrix.R"))
source(file.path(scripts.dir, funcs.dir, "MergeToLong.R"))
source(file.path(scripts.dir, funcs.dir, "PlotGeneAcrossTissues.R"))
source(file.path(scripts.dir, funcs.dir, "ConvertLongToWide.R"))
source(file.path(scripts.dir, funcs.dir, "ProjectToFrequency.R"))
source(file.path(scripts.dir, funcs.dir, "SvdFunctions.R"))  # many script-specific functions here

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
normalized.array.fname <- "array.adj.0.07.txt"
normalized.array.path <- file.path(data.dir, normalized.array.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)

# Load file ---------------------------------------------------------------

normalized.array <- LoadNormalizedArray(normalized.array.path)
rna.seq.exprs <- LoadRnaSeq(rna.seq.path)


# Log2 transform of array and rnaseq --------------------------------------

normalized.array <- log2(normalized.array + 1)
rna.seq.exprs <- log2(rna.seq.exprs + 1)


# Remove genes from RNA-Seqs that are not expressing. ---------------------

# # plot and show cutoff visually
# plot(density(unlist(rna.seq.exprs)))
# abline(v = cutoff, col = 'red')
# 
# jmean <- apply(rna.seq.exprs, 1, mean)
# 
# genes.not.expressed <- names(jmean[which(jmean < cutoff)])
# all.genes <- rownames(rna.seq.exprs)
# genes.expressed <- all.genes[which(!all.genes %in% genes.not.expressed)]
# 
# print(paste0(length(genes.not.expressed), 
#              " not expressed at cutoff = mean(", cutoff, 
#              ") for a gene across all tissues. Removing..."))
# 
# # Remove unexpressed genes
# 
# rna.seq.exprs.filtered <- rna.seq.exprs[genes.expressed, ]
rna.seq.exprs.filtered <- rna.seq.exprs  # no filter

# Take only common genes --------------------------------------------------

array.genes <- rownames(normalized.array)
rna.seq.genes <- rownames(rna.seq.exprs.filtered)

common.genes <- intersect(array.genes, rna.seq.genes)

print(paste(length(common.genes), "common genes between array and rnaseq"))

normalized.array <- normalized.array[common.genes, ]
rna.seq.exprs.filtered <- rna.seq.exprs.filtered[common.genes, ]

# Merge data into long format ---------------------------------------------

dat <- MergeToLong(normalized.array, rna.seq.exprs.filtered)

dat.split <- split(dat, dat$tissue)


# Remove WFAT from analysis -----------------------------------------------

# WFAT seems to be show a lot of epidydmal-specific genes, let's remove 

# it from analysis to prevent strangely oscillating genes showing up
dat.split$WFAT <- NULL

# remove BFAT from analysis: we see skeletal muscle and muscle contraction genes
# in PCA 2 (BFAT + Muscle)
# dat.split$BFAT <- NULL


# Project my data ---------------------------------------------------------

omega <- 2 * pi / 24
omegas <- GetOmegas()

start.time <- Sys.time()
if (getDoParWorkers() == 1){
  registerDoMC(40)
}
print(paste("Parallel processing on:", getDoParWorkers(), "cores"))
dat.split.proj <- lapply(dat.split, function(x){
  ddply(x, .(gene), ProjectToFrequency, my.omega = omega, normalize = FALSE, rhythmic.only = FALSE, pval.cutoff = 1, .parallel = TRUE)
})
print(Sys.time() - start.time)


# Add tissue information into each list -----------------------------------

for (tissue in names(dat.split.proj)){
  dat.split.proj[[tissue]]$tissue <- tissue
}

# Combine data ------------------------------------------------------------

dat.proj <- do.call(rbind, dat.split.proj)


# Remove NaNs -------------------------------------------------------------

# dat.proj <- dat.proj[is.finite(dat.proj$exprs.transformed), ]

# Clean up ----------------------------------------------------------------

# rm(dat.split.proj, dat.split)

# Do SVD on this normalized dataset ---------------------------------------

# long to wide conversion
dat.wide <- ConvertLongToWide(dat.proj, measurement.var = "exprs.transformed")

# Complete cases only. This removes NaN rows.
dat.wide <- dat.wide[complete.cases(dat.wide), ]


# Optional: centre rows ---------------------------------------------------

# # Centre the rows
# #dat.wide <- t(scale(t(dat.wide), center = TRUE, scale = FALSE))  # don't scale that's bad.
# 
# # Remove rows with avg 0
# dat.wide <- dat.wide[complete.cases(dat.wide), ]

s <- svd(dat.wide)


# Set rownames for u and v ------------------------------------------------

rownames(s$u) <- rownames(dat.wide)
rownames(s$v) <- colnames(dat.wide)

# Plot interesting components ---------------------------------------------

tissues <- GetTissues(colnames(normalized.array))

sing.vals <- seq(length(s$d))
# sing.vals <- c(1, 2, 3)

print("Plotting components...")

pdf("plots/components.WFAT.removed.pdf")

# ScreePlot
plot(s$d^2 / sum(s$d ^ 2), type='o')  # Manual screeplot. 

sink(file.path("results", "components.centred.WFAT.BFAT.removed.txt"))

for (sing.val in sing.vals){
  eigengene <- s$v[, sing.val]
  jmax <- max(Mod(eigengene))
  max.loading <- names(eigengene[which(Mod(eigengene) == jmax)])
  jmax.arg <- Arg(eigengene[max.loading])  # aligns everything relative to jmax.arg
  PlotComplex(eigengene, 
              axis.min = -jmax,
              axis.max = jmax,
              labels = names(eigengene), 
              col = "HSV",
              add.text.plot = FALSE, 
              main = paste("Component:", sing.val),
              jpch = 20,
              threshold = 0, 
              rotate = -jmax.arg)
  
  eigensample <- s$u[, sing.val]
  jmax <- max(Mod(eigensample))
  # rotate opposite way to make it all kosher
  PlotComplex(eigensample,
              axis.min = -jmax,
              axis.max = jmax,
              labels = names(eigensample),
              col = "HSV",
              main = paste("Component:", sing.val),
              add.text.plot = FALSE,
              jpch = 1,
              threshold = 0.5 * jmax,
              verbose = FALSE,
              rotate = jmax.arg)
  
  # show how fast the values drop
  top.genes.all <- GetTopNValues(Mod(eigensample), N = length(eigensample))
  plot(top.genes.all$vals, 
       log = "x", 
       main = paste("Component:", sing.val),
       ylab = "Modulus of gene loading",
       xlab = "Gene index")
  # Plot args in color
  top.genes <- GetTopNValues(Mod(eigensample), N = 100)# list of $vals $i
  # use drop to keep rownames
  outer.prod.mat <- s$d[sing.val] * OuterComplex(s$u[top.genes$i, sing.val, drop = FALSE], t(s$v[, sing.val, drop = FALSE]))
  outer.prod.mat <- OrderPhaseMatrix(outer.prod.mat, order.by = max.loading, order.tissues = TRUE)
  PlotArgsMatrix(outer.prod.mat, main = paste("Component", sing.val))
  
  # print top 100 genes (copy and paste-able) to file
  cat((paste0(c("COMPONENT", sing.val), collapse = "")))
  cat("\n")
  cat(paste0(head(names(top.genes$vals), n = 100), collapse = "\n"))
  cat("\n")
  
# #   # print top 20 genes (list-able)
#   cat(paste(head(names(top.genes$vals), n = 20), collapse = '", "'))
#   cat("\n")
}
sink()
dev.off()

# plot rhythmic genes

print("Plotting gene loadings")

pdf("plots/rhythmic.genes.by.components.WFAT.removed.pdf")
for (sing.val in sing.vals){
  eigensample <- s$u[, sing.val]
  jtissues <- rownames(s$v)
  top.genes <- GetTopNValues(Mod(eigensample), N = 10)# list of $vals $i
  
  genes <- names(top.genes$vals)
  
  dat.sub <- subset(dat, gene %in% genes)
  dat.sub$gene <- factor(dat.sub$gene)
  
  for (jgene in names(top.genes$vals)){
    print(PlotGeneAcrossTissues(subset(dat.sub, gene == jgene & tissue %in% jtissues), jtitle = paste("gene:", jgene, "component", sing.val)))
  }
}
dev.off()
