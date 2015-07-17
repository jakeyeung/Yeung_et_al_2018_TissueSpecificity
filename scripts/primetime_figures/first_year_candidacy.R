# 2015-07-10
# Jake Yeung
# Figures for first year candidacy

# library(Sushi)
setwd("~/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots"

# Functions ---------------------------------------------------------------

library(hash)
library(dplyr)
library(ggplot2)
library(grid)
library(gplots)
source("scripts/functions/LoadKallisto.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/TissueSpecificRhythmicsFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("~/projects/tissue-specificity/scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/DifferentialSitecountsFunctions.R")
source("scripts/functions/LoadSitecounts.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LongToMat.R")
source("scripts/functions/ShannonEntropy.R")
source("scripts/functions/FitRhythmic.R")

GetActSvd <- function(act.long, pval.adj.cutoff = 0.05){
  # Which ones are most rhythmic? -------------------------------------------
  act.fit <- act.long %>%
    group_by(tissue, gene) %>%
    do(FitRhythmicWeighted(dat = .))
  
  # False discovery rate adj ------------------------------------------------
  
  act.fit$pval.adj <- p.adjust(act.fit$pval, method = "BH")
  
  # Show top genes for each tissue ------------------------------------------
  
  head(arrange(data.frame(act.fit), pval), n = 50)
  
  library(plyr)
  source("~/projects/tissue-specificity/scripts/functions/SvdFunctions.R")  # many script-specific functions here
  
  omega <- 2 * pi / 24
  
  start.time <- Sys.time()
  act.complex <- lapply(split(act.long, act.long$tissue), function(x){
    ddply(x, .(gene), ProjectToFrequency2, omega = omega, add.tissue = TRUE)
  }) %>%
    do.call(rbind, .) %>%
    mutate(magnitude = Mod(exprs.transformed)) %>%
    arrange(desc(magnitude))
  print(Sys.time() - start.time)
  # head(act.complex, n = 100)
  
  detach("package:plyr", unload=TRUE)
  
  library(dplyr)
  
  # SVD on complex matrix ---------------------------------------------------
  # optinally filter out motifs
  
  filter.motifs <- subset(act.fit, pval.adj <= pval.adj.cutoff)$gene
  # act.complex.mat <- dcast(data = act.complex, formula = gene ~ tissue, value.var = "exprs.transformed")
  act.complex.mat <- dcast(data = subset(act.complex, gene %in% filter.motifs), formula = gene ~ tissue, value.var = "exprs.transformed")
  rownames(act.complex.mat) <- act.complex.mat[, "gene"]
  act.complex.mat <- act.complex.mat[, 2:ncol(act.complex.mat)]
  
  act.svd <- svd(act.complex.mat) 
  
  # add row and colnames
  rownames(act.svd$u) <- rownames(act.complex.mat)
  rownames(act.svd$v) <- colnames(act.complex.mat)
  return(act.svd)
}

CorrelatePromoterUsageToAmp <- function(tpm.filt, dat.rhyth.relamp, avgexprsdic, filter.tissues, tissue.spec){
  source("scripts/functions/AlternativeFirstExonsFunctions.R")
  
#   tpm.filt <- subset(tpm.filt, !tissue %in% filter.tissues)
  # tpm.filt <- tpm.filt
  tpm.afe <- tpm.filt %>%
    group_by(gene_name, tissue, time) %>%
    mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1))
  
  tested.genes <- unique(dat.rhyth.relamp$gene)  # 20013 genes
  
  tpm.avg <- tpm.afe %>%
    subset(., !tissue %in% filter.tissues & gene_name %in% tested.genes) %>%
    group_by(gene_name, transcript_id, tissue) %>%
    summarise(tpm_norm.avg = mean(tpm_normalized), tpm_norm.var = var(tpm_normalized))
  
  tpm.avg$relamp <- mapply(function(tiss, gene) relampdic[[as.character(paste(tiss, gene, sep=";"))]],
                           as.character(tpm.avg$tissue), 
                           as.character(tpm.avg$gene_name))
  
  tpm.avg$int.rnaseq <- mapply(function(tiss, gene) avgexprsdic[[as.character(paste(tiss, gene, sep=";"))]],
                               as.character(tpm.avg$tissue),
                               as.character(tpm.avg$gene_name))
  
  tpm.avg.filt <- subset(tpm.avg, gene_name %in% tissue.spec & tpm_norm.var > 0 & int.rnaseq > 4)
  
  tpm.fit <- tpm.avg.filt %>%
    group_by(gene_name, transcript_id) %>%
    do(FitPromoterUsageToAmplitude(.))
  
  tpm.summary <- tpm.fit %>%
    group_by(gene_name) %>%
    do(SubsetMinPval(jdf = .))
  
  tpm.summary$pval.adj <- p.adjust(tpm.summary$pval)  
  return(list(tpm.summary = tpm.summary, tpm.avg = tpm.avg, tpm.fit = tpm.fit))
}

figcount <- 1

# Figure 1: tissue-specific rhythms ---------------------------------------
# first show MARA results

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)
figcount <- figcount + 1

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long <- LoadActivitiesLong(indir)

act.svd <- GetActSvd(act.long, pval.adj.cutoff = 0.0005)

eigens <- GetEigens(act.svd, comp = 1)

# now show some tissue-specific results
# code from find_oscillating_genes.pairs.R
load("Robjs/dat.rhyth.relamp.pvalmin1e-5.pvalmax0.05.relampmax.0.1.meancutoff6.Robj", verbose = T)
load("Robjs/arrayrnaseq.mydat.Robj", verbose = T)

# how many rhythmic genes?
n.rhyth <- dat.rhyth.relamp %>%
  group_by(tissue) %>%
  summarise(n.rhyth = length(is.rhythmic[which(is.rhythmic == TRUE)]))

n.pairs <- dat.rhyth.relamp %>%
  group_by(gene) %>%
  do(RhythTiss(.))

n.pairs.counts <- n.pairs %>%
  group_by(n.tiss) %>%
  summarise(count = length(n.tiss))

n.pairs.tissuecounts <- n.pairs %>%
  filter(n.tiss == 2) %>%
  group_by(tissues) %>%
  summarise(count = length(tissues)) %>%
  mutate(n.tiss = sapply(tissues, function(x) length(strsplit(x, split = ",")[[1]]))) %>%
  arrange(desc(count)) %>%
  filter(count > 10)
n.pairs.tissuecounts$tissues <- factor(n.pairs.tissuecounts$tissues, n.pairs.tissuecounts$tissues)

# global tissue-specific rhythmic
n.tissues.rhyth.plot <- ggplot(subset(n.pairs.counts, n.tiss >= 1), aes(x = n.tiss, y = count)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete() +
  ggtitle("Number of tissues in which gene oscillates") +
  ylab("Number of genes") +
  xlab("Number of tissues that is rhythmic for gene")

pairs.plot <- ggplot(subset(n.pairs.tissuecounts, n.tiss == 2), aes(x = tissues, y = count)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete() +
  ggtitle("Genes rhythmic in pairs of tissues") +
  ylab("Number of genes") +
  xlab("Number of tissues that is rhythmic for gene") +
  theme(axis.text.x=element_text(angle=90,vjust = 0))

jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
multiplot(eigens$u.plot, eigens$v.plot, n.tissues.rhyth.plot, pairs.plot, layout = jlayout)

dev.off()


# Figure 1 redo: tissue-specific and tissue-wide rhythms ------------------
# from heatmap_tissue_specific_rhythms.R

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, "B_maxamp.pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)

load(file = "Robjs/dat.fit.Robj", verbose = TRUE)
load(file = "Robjs/dat.complex.maxexprs4.maxrelamp01.minpval1e-3.Robj", verbose = TRUE)

ref.gene <- "Arntl"
filt.tiss.lst <- list(c(), c("Cere", "BS", "Hypo"))

for (filt.tiss in filt.tiss.lst){
  dat.fit.relamp <- GetRelamp(dat.fit, max.pval = 1e-3)
  print(subset(dat.fit.relamp, relamp == 1))
  #   dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)

  
  dat.fit.relamp <- dat.fit.relamp %>%
    group_by(gene) %>%
    mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))
  
  dat.entropy <- dat.fit.relamp %>%
    #   subset(., min.pval <= 1e-5 & max.relamp > 0.1) %>%
    subset(., max.exprs >= 4 & max.relamp >= 0.1 & min.pval <= 1e-3 & !tissue %in% filt.tiss) %>%
    group_by(gene) %>%
    summarise(entropy = ShannonEntropy(relamp / sum(relamp)))
  
  N <- 2000
  n.genes <- nrow(dat.entropy)
  low.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy), ]$gene, n = N)
  high.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]$gene, n = N)
  
  M.low <- LongToMat(subset(dat.fit.relamp, gene %in% low.entropy.genes & !tissue %in% filt.tiss))
  M.high <- LongToMat(subset(dat.fit.relamp, gene %in% high.entropy.genes & !tissue %in% filt.tiss))
  
  dat.complex$real <- Re(dat.complex$exprs.transformed)
  dat.complex$imag <- Im(dat.complex$exprs.transformed)
  
  M.real <- t(LongToMat(subset(dat.complex, !tissue %in% filt.tiss), value.var = "real"))
  M.imag <- t(LongToMat(subset(dat.complex, !tissue %in% filt.tiss), value.var = "imag"))
  M.all <- cbind(M.real, M.imag)
  hc <- hclust(dist(M.all, method = "euclidean"))
  
  dat.ref <- subset(dat.fit.relamp, relamp == 1)
  dic <- hash(as.character(dat.ref$tissue), dat.ref$amp / 2)
  
  dat.complex <- dat.complex %>%
    group_by(tissue) %>%
    do(NormalizeComplexMat(., dic))
  
  dat.complex$real.norm <- Re(dat.complex$exprs.norm)
  dat.complex$imag.norm <- Im(dat.complex$exprs.norm)
  
  M.real.norm <- t(LongToMat(subset(dat.complex, !tissue %in% filt.tiss), value.var = "real.norm"))
  M.imag.norm <- t(LongToMat(subset(dat.complex, !tissue %in% filt.tiss), value.var = "imag.norm"))
  M.all.norm <- cbind(M.real.norm, M.imag.norm)
  hc.norm <- hclust(dist(M.all.norm, method = "euclidean"))
  
  par(mfrow=c(2,2))
#   plot(hc, main = "Clustering on fourier component T=24h")
  plot(density(dat.entropy$entropy[which(!is.na(dat.entropy$entropy))]), main = "Distribution of entropy in rhythmic genes")
  plot(hc.norm, main = "Clustering: fourier T=24h (normalized by max amp)")
  PlotRelampHeatmap(M.low, paste0("Low entropy rhythmic genes. N=", length(low.entropy.genes)))
  PlotRelampHeatmap(M.high, paste0("High entropy rhythmic genes. N=", length(high.entropy.genes)))
  
}
dev.off()

# Figure 2: kidney-liver and bfat-liver -----------------------------------


pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)
figcount <- figcount + 1

# BEGIN: Kidney to Liver comparison
genes.kidliv <- subset(n.pairs, tissues == "Kidney,Liver")
phases.kidliv <- subset(dat.rhyth.relamp, gene %in% genes.kidliv$gene & tissue %in% c("Kidney", "Liver"))
phases.kidliv.diff <- phases.kidliv %>%
  group_by(gene) %>%
  summarise(phasediff = GetPhaseDiff(phase))
kidney.liver.phasediff.plot <- ggplot(phases.kidliv.diff, aes(x = phasediff)) + scale_y_discrete() +
  geom_histogram(binwidth = 1) + ggtitle("Liver phase - Kidney phase")  # liver minus kidney shows coherence
kidney.liver.phase.plot <- ggplot(subset(phases.kidliv), aes(x = phase, fill = tissue)) + scale_y_discrete() +
  geom_histogram(binwidth = 1) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
# END: Kidney to Liver comparison

# BEGIN: Liver to BFAT comparison
genes.bfatliv <- subset(n.pairs, tissues == "BFAT,Liver")
phases.bfatliv <- subset(dat.rhyth.relamp, gene %in% genes.bfatliv$gene & tissue %in% c("BFAT", "Liver"))
phases.bfatliv.diff <- phases.bfatliv %>%
  group_by(gene) %>%
  summarise(phasediff = GetPhaseDiff(phase))
bfat.liver.phasediff.plot <- ggplot(phases.bfatliv.diff, aes(x = phasediff)) + scale_y_discrete() +
  geom_histogram(binwidth = 1) + ggtitle("Liver phase - BFAT phase")  # liver minus kidney shows coherence
bfat.liver.phase.plot <- ggplot(subset(phases.bfatliv), aes(x = phase, fill = tissue)) + scale_y_discrete() +
  geom_histogram(binwidth = 1) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
# END: Liver to BFAT comparison

# BEGIN: BFAT to Mus comparison
genes.bfatmus <- subset(n.pairs, tissues == "BFAT,Mus")
phases.bfatmus <- subset(dat.rhyth.relamp, gene %in% genes.bfatmus$gene & tissue %in% c("BFAT", "Mus"))
phases.bfatmus.diff <- phases.bfatmus %>%
  group_by(gene) %>%
  summarise(phasediff = GetPhaseDiff(phase))
bfat.mus.phasediff.plot <- ggplot(phases.bfatmus.diff, aes(x = phasediff)) + scale_y_discrete() +
  geom_histogram(binwidth = 1) + ggtitle("Mus phase - BFAT phase")  # liver minus kidney shows coherence
bfat.mus.phase.plot <- ggplot(subset(phases.bfatmus), aes(x = phase, fill = tissue)) + scale_y_discrete() +
  geom_histogram(binwidth = 1) + coord_polar(theta = "x") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
# END: Liver to BFAT comparison

jlayout2 <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
multiplot(kidney.liver.phasediff.plot, kidney.liver.phase.plot, bfat.liver.phasediff.plot, bfat.liver.phase.plot, layout = jlayout2)

dev.off()

# Figure 3: Hnf and Mef2 --------------------------------------------------

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)
figcount <- figcount + 1

hnf4a.activity <- ggplot(subset(act.long, gene == "HNF4A_NR2F1.2.p2" & experiment == "rnaseq"), 
                         aes(x = time, y = exprs)) +
  geom_line() +
  geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
  facet_wrap(~tissue) + 
  xlab("CT") +
  ylab("Activity") + 
  ggtitle("Hnf4a motif activity") +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))  +
  theme(axis.text.x=element_text(angle=90,vjust = 0))

hnf4a.exprs <- ggplot(subset(mydat, gene == "Hnf4a" & experiment == "rnaseq"), aes(x = time, y = exprs)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Hnf4a mRNA expression") + 
  ylab(label = "log2 mRNA expression") +
  xlab("CT") +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0))

hnf4a.regulated.examp <- ggplot(subset(mydat, gene == "Upp2" & experiment == "rnaseq"), aes(x = time, y = exprs)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Upp2 (regulated by Hnf4a)") + 
  ylab(label = "log2 mRNA expression") +
  xlab("CT") +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))  +
  theme(axis.text.x=element_text(angle=90,vjust = 0))

mef2c.activity <- ggplot(subset(act.long, gene == "MEF2.A.B.C.D..p2" & experiment == "rnaseq"), aes(x = time, y = exprs)) +
  geom_line() +
  geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
  facet_wrap(~tissue) + 
  xlab("CT") +
  ylab("Activity") + 
  ggtitle("Mef2 motif activity") +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0))

mef2c.exprs <- ggplot(subset(mydat, gene == "Mef2c" & experiment == "rnaseq"), aes(x = time, y = exprs)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Mef2c mRNA expression") + 
  ylab(label = "log2 mRNA expression") +
  xlab("CT")  +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0))

mef2c.regulated.examp <- ggplot(subset(mydat, gene == "Des" & experiment == "rnaseq"), aes(x = time, y = exprs)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Des (regulated by Mef2c)") + 
  ylab(label = "log2 mRNA expression") +
  xlab("CT")  +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0))

jlayout3 <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE)
multiplot(hnf4a.activity, hnf4a.exprs, hnf4a.regulated.examp, mef2c.activity, mef2c.exprs, mef2c.regulated.examp, layout = jlayout3)

dev.off()

# Figure 4: Alternative promoter ------------------------------------------

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)
figcount <- figcount + 1

load(file = "Robjs/tpm.merged.Robj", verbose = T)

# for annotating rhythmic genes
keys <- paste(dat.rhyth.relamp$tissue, dat.rhyth.relamp$gene, sep = ";")
relampdic <- hash(keys, dat.rhyth.relamp$relamp)
avgexprsdic <- hash(keys, dat.rhyth.relamp$int.rnaseq)

tissue.spec <- unique(subset(dat.rhyth.relamp, is.tiss.spec == TRUE)$gene)
filter.tissues <- c()
tpm.filt <- subset(tpm.merged, !tissue %in% filter.tissues & gene_name %in% tissue.spec)

# ~50 seconds
start.time <- Sys.time()
tpm.output <- CorrelatePromoterUsageToAmp(tpm.filt, dat.rhyth.relamp, avgexprsdic, filter.tissues, tissue.spec)
tpm.summary <- tpm.output$tpm.summary
tpm.avg <- tpm.output$tpm.avg
tpm.fit <- tpm.output$tpm.fit
print(Sys.time() - start.time)

pval.cutoff <- 0.05
tpm.summary.filt <- subset(tpm.summary, pval <= pval.cutoff)
sig.hits <- tpm.summary.filt$gene_name

# plot examples or top hits
jgene <- "Upp2"
tpm.sub <- subset(tpm.fit, gene_name == jgene)
jtranscript <- tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1]
upp2.altprom.plot <- ggplot(subset(tpm.avg, transcript_id == jtranscript), aes(y = tpm_norm.avg, x = relamp, label = tissue)) + 
  geom_point() + geom_text() + ggtitle(jgene) + geom_smooth(method = "lm") + xlab("Relative amplitude") + ylab("Fraction promoter usage")

# plot histogram of p-values
altprom.histo <- ggplot(tpm.fit, aes(x = pval)) + geom_histogram(binwidth = 0.0075) + xlab("P-value")

jlayout4 <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
multiplot(altprom.histo, upp2.altprom.plot, layout = jlayout4)

dev.off()


# Figure 5: DHS and rhythmicity -------------------------------------------

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)
figcount <- figcount + 1

N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_sum_matrix_bugfixed"
suffix <- "dist_filt.bugfixed.sitecount.matrix"
N <- LoadSitecountsEncodeAll(maindir = N.dir, tissues = c("Liver", "Kidney", "Cere", "Lung", "Heart", "Mus"),
                             suffix = suffix, with.ensemblid = FALSE, rename.tissues = FALSE)  # merged by gene
dhs.tiss <- unique(N$tissue)
tfs <- GetTFs()

jgene <- c("S100a10")

PlotBiplot <- function(N, jgene){
  X.motif <- GetMotifMat(N, jgene = jgene)
  X.motif <- t(scale(t(X.motif), center = TRUE, scale = FALSE))
  s <- prcomp(X.motif, center = FALSE, scale. = FALSE)
  bplot <- PCbiplot(PC = s, jtitle = paste(jgene, "DHS site counts"))
  return(bplot) 
}

# jgenes <- c("S100a10", "Fmo5", "Amt", "Sigmar1", "Mcm4")
jgenes <- c("Ddc", "Pnp", "Insig2")

exprs.plots <- lapply(jgenes, function(g) return(PlotRnaseqAcrossTissues(subset(mydat, gene == g & experiment == "rnaseq"), jtitle = paste(g, "RNA-seq profile"))))
biplots <- lapply(jgenes, function(g) return(PlotBiplot(N, g)))

jlayout5 <- matrix(c(1, 4, 2, 5, 3, 6), 3, 2, byrow = TRUE)
do.call(multiplot, c(exprs.plots, biplots, list(layout = jlayout5)))

dev.off()