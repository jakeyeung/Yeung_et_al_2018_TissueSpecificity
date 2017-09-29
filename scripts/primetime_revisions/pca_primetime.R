# Jake Yeung
# Date of Creation: 2017-09-22
# File: ~/projects/tissue-specificity/scripts/primetime_figures/pca_primetime.R
# Plot PCAs in primetime style. In 2D, over time, with gene loadings to illustrate temporal and tissue variation (maybe show a global mean)

rm(list=ls())

remove.wfat <- TRUE

library(ggplot2)
library(hash)
library(dplyr)

# Define my functions -----------------------------------------------------

# First source my functions I wrote
funcs.dir <- file.path('scripts', 'functions')
source(file.path(funcs.dir, 'SampleNameHandler.R'))  # for shortening sample names
source(file.path(funcs.dir, 'PcaPlotFunctions.R'))  # for visualizing PCA, periodoigrams
source(file.path(funcs.dir, 'FourierFunctions.R'))  # for periodoigrams
source(file.path(funcs.dir, 'GetTissueSpecificMatrix.R'))  # as name says
source(file.path(funcs.dir, "GrepRikGenes.R"))
source(file.path(funcs.dir, "GetTissueTimes.R"))
source(file.path(funcs.dir, "PCAFunctions.R"))
source(file.path(funcs.dir, "LoadAndHandleData.R"))
source(file.path(funcs.dir, "FitRhythmic.R"))
source(file.path(funcs.dir, "PlotGeneAcrossTissues.R"))
source(file.path(funcs.dir, "LoadArray.R"))

BinVector <- function(x, ...){
  # bin vector into discrete bins i
  ends <- c(...)
  starts <- c(1, ends[1:length(ends) - 1])
  bins <- sapply(x, function(y){
    for (i in seq(length(starts))){
      if (y >= starts[i] & y <= ends[i]){
        return(i)
      } else {
        next
      }
    }
  })
  return(unlist(bins))
}

library(scales)
mylog_trans <- function(base=exp(1), from=0){
  # from
  # http://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), 
            domain = c(base^from, Inf))
}

PlotGeneLoadings <- function(x, jpc1, jpc2, top.n){
  genesx <- rownames(x)[order(x[, jpc1] ^ 2, decreasing = TRUE)][1:top.n]
  genesy <- rownames(x)[order(x[, jpc2] ^ 2, decreasing = TRUE)][1:top.n]
  gene.hash <- hash(c(genesx, genesy), TRUE)
  genes.lab <- sapply(rownames(x), function(g) ifelse(!is.null(gene.hash[[g]]), g, NA))
  jdf <- data.frame(pcx = x[, jpc1], pcy = x[, jpc2], gene = genes.lab)
  m <- ggplot(jdf, aes(x = pcx, y = pcy, label = gene)) + 
    geom_point(alpha = 0.1) + 
    xlab(jpc1) + ylab(jpc2) + 
    geom_text_repel(size = 7) + 
    theme_bw() + 
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  return(m)
}

PlotGeneLoadings.1d <- function(x, jpc1, top.n){
  genesx <- order(x[, jpc1] ^ 2, decreasing = TRUE)[1:top.n]
  x.sub <- sort(x[genesx, jpc1])
  jdf <- data.frame(pcx = x.sub, gene = names(x.sub), stringsAsFactors=FALSE)
  jdf$gene <- factor(jdf$gene, levels = names(x.sub))
  m <- ggplot(jdf, aes(x = gene, y = pcx)) + 
    geom_bar(stat = "identity") +
    xlab(jpc1) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
}


PlotSampleLoadings <- function(pca.long, jpc, draw.global.mean = FALSE, draw.tissue.mean = FALSE, tissue.order = "auto"){
  jsub <- subset(pca.long, pc == jpc)
  
  if (tissue.order == "auto"){
    jsub.mean <- jsub %>%
      group_by(tissue) %>%
      summarise(loading = mean(loading)) %>%
      arrange(desc(loading))
    tissue.order <- as.character(jsub.mean$tissue)
  }
  if (!is.null(tissue.order)){
    jsub$tissue <- factor(as.character(jsub$tissue), levels = tissue.order)
  }
  
  m <- ggplot(jsub, aes(x = time, y = loading)) + 
    geom_point() + geom_line() + facet_wrap(~tissue) + 
    theme_bw() + 
    ggtitle(paste("Principal Component:", jpc)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  if (draw.global.mean){
    g.mean <- mean(jsub$loading)
    m <- m + geom_hline(yintercept = g.mean, linetype = "dotted")
  }
  if (draw.tissue.mean){
    jsub.means <- jsub %>%
      group_by(tissue) %>%
      summarise(jmean = mean(loading))
    m <- m + geom_hline(data = jsub.means, mapping = aes(yintercept = jmean), linetype = "dotted")
  }
  m <- m + xlab("Time (CT)") + ylab("Loadings [A.U.]") + 
    scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))
  return(m)
}

PlotLoadingsOneTissue <- function(pca.long, jtiss, pcs){
  pca.sub <- subset(pca.long, tissue == jtiss & pc %in% pcs)
  m <- ggplot(pca.sub, aes(x = time, y = loading)) + 
    geom_point() + geom_line() + 
    theme_bw() + 
    ggtitle(paste0(jtiss, " loadings")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~pc) + 
    xlab("Time (CT)") + ylab("Loadings [A.U.]") + 
    scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) + 
    geom_hline(yintercept = 0, linetype = "dotted")
  return(m)
}

# Set constants -----------------------------------------------------------

outdir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots_full_paper_revised"
dir.create(outdir)

# Figure 1A Present high level variation BETWEEN TISSUES and WITHIN TISSUES (time) ------------------
# scripts from pca_adjusted_microarray.label_variance.R

# Load data, log transform

log2.transform <- TRUE
# load("Robjs/dat.array.adj.primetime.Robj", verbose=T)  # no -Infs

dat <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt",
                           remove.negs = TRUE, fix.rik.xgene = TRUE)

# Transform

if (log2.transform){
  dat <- log(dat + 1, base = 2)
}

# no WFAT
if (remove.wfat){
  dat <- dat[, !grepl("WFAT", colnames(dat))] 
}

# PCA center by ROW
dat.centered <- dat - apply(dat, 1, mean)
dat_pca <- prcomp(t(dat.centered), center=FALSE, scale.=FALSE)
sdev.norm <- sapply(dat_pca$sdev, function(x) x ^ 2 / sum(dat_pca$sdev ^ 2))


# Consolidate PCA into long
jtissues <- GetTissues(rownames(dat_pca$x), get_unique = FALSE)
jtimes <- as.numeric(GetTimes(rownames(dat_pca$x), get_unique = FALSE))
pca.long <- data.frame(tissue = rep(jtissues, times = ncol(dat_pca$x)),
                       time = rep(jtimes, times = ncol(dat_pca$x)),
                       pc = rep(colnames(dat_pca$x), each = nrow(dat_pca$x)),
                       loading = unlist(as.vector(dat_pca$x), use.names = FALSE))
head(pca.long)
n.samps <- length(jtissues)

plot.i <- 0
pdf(file.path(outdir, paste0(plot.i, ".PCA_component_vs_component.pdf")), useDingbats=FALSE)

# Plot PC1 vs PC2
jpc1 <- "PC1"
jpc2 <- "PC2"
x <- subset(pca.long, pc == jpc1)$loading
y <- subset(pca.long, pc == jpc2)$loading
tisslab <- subset(pca.long, pc == jpc1)$tissue
timelab <- subset(pca.long, pc == jpc1)$time
cols.uniq <- rainbow(length(unique(tisslab)))
jpch <- as.numeric(seq(length(unique(tisslab))))
jpch[which(jpch == 11)] <- 20
jpch[4] <- 44
jpch[5] <- 96
jpch[7] <- 46
jpch.vec <- rep(jpch, each = length(unique(timelab)))

cols <- rep(cols.uniq, each = length(unique(timelab)))
tisstimelab <- paste(as.character(tisslab), as.character(timelab))
# with labels
# plot(x, y, cex = 0.01, xlab = jpc1, ylab = jpc2)
# text(x, y, labels = tisstimelab, col = cols)  
# no labels
plot(x, y, cex = 1.7, main = paste0(jpc1, " vs. ", jpc2), xlab = jpc1, ylab = jpc2, pch = jpch.vec, bg = "white")
# legend("bottomright", as.character(unique(tisslab)), pch = 19, title = "Tissue", col = cols.uniq, horiz = F)
legend("bottomright", as.character(unique(tisslab)), title = "Tissue", pch = jpch, horiz = F)

# Plot PC13 vs PC17
library(PhaseHSV)
jtiss <- "Liver"
jpc1 <- "PC13"
jpc2 <- "PC17"
x <- subset(pca.long, pc == jpc1 & tissue == jtiss)$loading
y <- subset(pca.long, pc == jpc2 & tissue == jtiss)$loading
time <- as.numeric(subset(pca.long, pc == jpc1 & tissue == jtiss)$time)
time.mod <- time %% 24
time.cols <- hsv(PhaseToHsv(2 * pi * time.mod / 24, 0, 2 *pi), s=0.9, v=0.7)
cols <- as.numeric(subset(pca.long, pc == jpc1)$tissue)
plot(x, y, type = "n", xlab = jpc1, ylab = jpc2)
text(x, y, labels = time, col = time.cols, cex = 2)  

dev.off()


# Plot screeplot with tissue/temporal variance ----------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch
dat.long <- subset(dat.long, tissue != "WFAT")

top.n <- 10

source("/home/yeung/projects/tissue-specificity/scripts/functions/ListFunctions.R")

plst <- expandingList()
# tissue variance
jpc1 <- "PC1"; jpc2 <- "PC2"
plst$add(PlotSampleLoadings(pca.long, jpc1, draw.global.mean = TRUE, tissue.order = "auto"))
plst$add(PlotSampleLoadings(pca.long, jpc2, draw.global.mean = TRUE, tissue.order = "auto"))
plst$add(PlotGeneLoadings(dat_pca$rotation, jpc1, jpc2, top.n))
plst$add(PlotGeneAcrossTissues(subset(dat.long, gene == "Gabra1" & experiment == "array"), make.pretty = TRUE))
plst$add(PlotGeneAcrossTissues(subset(dat.long, gene == "Cox8b" & experiment == "array"), make.pretty = TRUE))

# temporal variance
jpc1 <- "PC13"; jpc2 <- "PC17"
plst$add(PlotSampleLoadings(pca.long, jpc1, draw.tissue.mean = TRUE, tissue.order = "auto"))
plst$add(PlotSampleLoadings(pca.long, jpc2, draw.tissue.mean = TRUE, tissue.order = "auto"))
plst$add(PlotGeneLoadings(dat_pca$rotation, jpc1, jpc2, top.n))
plst$add(PlotGeneAcrossTissues(subset(dat.long, gene == "Arntl" & experiment == "array"), make.pretty = TRUE))
plst$add(PlotGeneAcrossTissues(subset(dat.long, gene == "Serpina1c" & experiment == "array"), make.pretty = TRUE))

# temporal variance PC13 to PC20
pcs <- paste("PC", seq(13, 21), sep = "")
jtiss <- "Liver"
plst$add(PlotLoadingsOneTissue(pca.long, jtiss, pcs))
jtiss <- "Kidney"
plst$add(PlotLoadingsOneTissue(pca.long, jtiss, pcs))

plst <- plst$as.list()
pdf(file.path(outdir, paste0(plot.i, ".PCA_loadings_gene_and_samps.pdf")), useDingbats=FALSE)
  lapply(plst, print)
dev.off()


# Analyze whether each PC is rhythmic or not ------------------------------

pca.p <- pca.long %>%
  group_by(pc, tissue) %>%
  do(GetPeriodogramFreq(.))

# Summarize by number of tissues with 12, 24, Inf, or other ---------------

pca.p.sum <- pca.p %>%
  group_by(pc) %>%
  do(SummarisePeriodogram2(., weighted = TRUE))

pca.p.sum$pc.num <- sapply(pca.p.sum$pc, function(p) as.numeric(strsplit(as.character(p), split = "PC")[[1]][[2]]))
pca.p.sum <- pca.p.sum[order(pca.p.sum$pc.num), ]
head(data.frame(pca.p.sum), n = 25)

# Plot screeplot, label each component by number of tissues in eac --------

eigenvals <- dat_pca$sdev ^ 2 / sum(dat_pca$sdev ^ 2)
pcs <- c(paste("PC", seq(n.samps), sep = ""))
eigenvals.dic <- hash(pcs, eigenvals)

# adjust factors for plotting
pca.p.sum$eigenvals <- sapply(as.character(pca.p.sum$pc), function(pc) eigenvals.dic[[pc]])
pca.p.sum$pc.num <- sapply(as.character(pca.p.sum$pc), function(x) as.numeric(substr(x, 3, nchar(x))))
pca.p.sum$pc <- factor(as.character(pca.p.sum$pc), levels = pcs)

# Make long and plot -----------------------------------------------------

pca.p.sum.long <- melt(pca.p.sum, id.var = c("pc", "eigenvals", "pc.num"), value.name = "fracFourier", variable.name = "Component")

pca.p.sum.long <- pca.p.sum.long %>%
  mutate(eigenvals.frac = eigenvals * fracFourier)

pca.p.sum.long$Component <- factor(pca.p.sum.long$Component, levels = c("T.inf", "T.other",  "T.12", "T.24"))
pca.p.sum.long <- pca.p.sum.long %>%
  group_by(pc) %>%
  arrange(Component)

starts <- c(1, 13, 100)
ends <- c(12, 99, 288)

pdf(file.path(outdir, paste0(plot.i, ".PCA_eigenvalues.pdf")), useDingbats=FALSE)
for (i in seq(length(starts))){
  print(PlotComponents2(pca.p.sum.long, starts[[i]], ends[[i]]))
}
dev.off()


# # Parseval's Theorem (check that our calculations are correct) ------------------------------------------------------
# 
# pca.p.parseval <- pca.p %>%
#   group_by(pc) %>%
#   summarise(eigen.est = sum(p.inf + p.time) / length(p.inf))
# pca.p.parseval$pc.num <- sapply(pca.p.parseval$pc, function(p) as.numeric(strsplit(as.character(p), "PC")[[1]][[2]]))
# pca.p.parseval <- pca.p.parseval[order(pca.p.parseval$pc.num), ]
# 
# eigenvals <- dat_pca$sdev ^ 2
# 
# diff.norm <- abs(eigenvals - pca.p.parseval$eigen.est) ^ 2 / eigenvals
# 
# plot(density(diff.norm))
# 
# plot(pca.p.parseval$eigen.est, eigenvals)  # should be close!
# 
# # Plot all PCs across conditions with gene loadings -----------------------
# 
# # show PC1 and PC2 gene loadings
# 
# 
# # get gene loadings into long format
# jpc <- "PC13"
# sorted.eigensamp <- dat_pca$rotation[, jpc][order(abs(dat_pca$rotation[, jpc]), decreasing = TRUE)]
# sorted.eigengene <- dat_pca$x[, jpc][order(abs(dat_pca$x[, jpc]), decreasing = TRUE)]
# head(sorted.eigengene)
# plot(dat_pca$x[, jpc])
# maxi <- 30
# 
# pcs <- seq(50); pcs <- paste("PC", pcs, sep = "")
# 
# pdf(file.path(outdir, paste0(plot.i, ".PCA_tissue_over_time.pdf")), useDingbats=FALSE)
# for (jpc in pcs){
#   print(ggplot(subset(pca.long, pc == jpc), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jpc))
# }
# dev.off()

