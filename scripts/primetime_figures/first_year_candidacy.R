# 2015-07-10
# Jake Yeung
# Figures for first year candidacy

# library(Sushi)
setwd("~/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots2"
dir.create(outdir)

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
source("scripts/functions/GetActSvd.R")
source("scripts/functions/SortByTissue.R")

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

# Figure 1 redo: SVD analysis ---------------------------------------------

load(file = "Robjs/dat.complex.maxexprs4.Robj")
load(file = "Robjs/dat.fit.Robj")
load(file = "Robjs/dat.long.Robj")

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, "C.pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)


ref.gene <- "Nr1d1"

# dat.fit.relamp <- GetRelamp(dat.fit, max.pval = 1e-3)
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)

dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))

# dat.complex <- TemporalToFrequencyDatLong(subset(dat.long, gene %in% dat.entropy$gene), period = 24, n = 8, interval = 6, add.entropy.method = "array")
# # adjust for "noise" and normalize to reference gene (Nr1d1)
# ref.amps <- subset(dat.complex, gene == ref.gene, select = c(tissue, gene, mag.norm))
# ref.amps.dic <- hash(ref.amps$tissue, ref.amps$mag.norm)
# dat.complex$exprs.adj <- dat.complex$exprs.transformed * dat.complex$frac.weight
# dat.complex$exprs.adj.norm <- mapply(function(tiss, exprs) exprs / ref.amps.dic[[tiss]], as.character(dat.complex$tissue), dat.complex$exprs.adj)
# dat.complex$mod.exprs.adj <- Mod(dat.complex$exprs.adj)
# dat.complex$mod.exprs.adj.norm <- Mod(dat.complex$exprs.adj.norm)

genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)

# ref gene should be max entropy as a check
dat.entropy <- dat.complex %>%
  subset(., gene %in% genes.exprs) %>%
  group_by(gene) %>%
  summarise(entropy = ShannonEntropy(Mod(exprs.adj.norm) / sum(Mod(exprs.adj.norm))))

# split genes into 3rds
N <- floor(nrow(dat.entropy) / 3)
n.genes <- nrow(dat.entropy)
low.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy), ]$gene, n = N)
high.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]$gene, n = N)
med.entropy.genes <- dat.entropy[order(dat.entropy$entropy), ]$gene[N:(n.genes - N)]

s.low.norm <- SvdOnComplex(subset(dat.complex, gene %in% low.entropy.genes), value.var = "exprs.adj")
s.high.norm <- SvdOnComplex(subset(dat.complex, gene %in% high.entropy.genes), value.var = "exprs.adj")

jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)

eigens.high.norm <- GetEigens(s.high.norm, period = 24, comp = 1)
eigens.high.examp <- PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Dbp" & experiment == "rnaseq"))
eigens.low1.norm <- GetEigens(s.low.norm, period = 24, comp = 1)
eigens.low1.examp <- PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Rgs16" & experiment == "rnaseq"))
eigens.low2.norm <- GetEigens(s.low.norm, period = 24, comp = 2)
eigens.low2.examp <- PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Myh7" & experiment == "rnaseq"))
eigens.low3.norm <- GetEigens(s.low.norm, period = 24, comp = 3)
eigens.low3.examp <- PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Ms4a1" & experiment == "rnaseq"))
# multiplot(eigens.high.norm$v.plot, eigens.high.norm$u.plot, eigens.high.examp, 
#           eigens.low1.norm$v.plot, eigens.low1.norm$u.plot, eigens.low1.examp, 
#           layout = jlayout)
multiplot(eigens.high.norm$v.plot, eigens.high.norm$u.plot,
          eigens.low1.norm$v.plot, eigens.low1.norm$u.plot, layout = jlayout)
multiplot(eigens.low2.norm$v.plot, eigens.low2.norm$u.plot, 
          eigens.low3.norm$v.plot, eigens.low3.norm$u.plot, layout = jlayout)
multiplot(eigens.high.examp, eigens.low1.examp, eigens.low2.examp, eigens.low3.examp, layout = jlayout)

dev.off()

# Figure 1 redo: tissue-specific and tissue-wide rhythms ------------------
# from heatmap_tissue_specific_rhythms.R

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, "B.pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)

filt.tiss <- c()

# dat.fit.relmap, dat.complex obtained above
M.low <- LongToMat(subset(dat.complex, gene %in% low.entropy.genes & !tissue %in% filt.tiss), value.var = "mod.exprs.adj.norm")
M.high <- LongToMat(subset(dat.complex, gene %in% high.entropy.genes & !tissue %in% filt.tiss), value.var = "mod.exprs.adj.norm")

dat.complex$real <- Re(dat.complex$exprs.adj)
dat.complex$imag <- Im(dat.complex$exprs.adj)

M.real <- t(LongToMat(subset(dat.complex, !tissue %in% filt.tiss), value.var = "real"))
M.imag <- t(LongToMat(subset(dat.complex, !tissue %in% filt.tiss), value.var = "imag"))
hc <- hclust(dist(cbind(M.real, M.imag), method = "euclidean"))

# summarize by spectral power
dat.mag <- dat.complex %>%
  group_by(tissue) %>%
  summarise(power = sum(mod.exprs.adj ^ 2)) %>%
  arrange(desc(power))
dat.mag$tissue <- factor(dat.mag$tissue, dat.mag$tissue)

# # Plot cluster
# M.args <- t(LongToMat(dat.complex, value.var = "exprs.transformed"))
# hc.args <- hclust(dist(Arg(M.args), method = "euclidean"))
# plot(hc.args)

# par(mfrow=c(2,2))
# #   plot(hc, main = "Clustering on fourier component T=24h")
# plot(density(dat.entropy$entropy[which(!is.na(dat.entropy$entropy))]), main = "Distribution of entropy in rhythmic genes")
# # draw 2 vertical lines for specifying tissue-spec to tissue-wide
# abline(v = sort(dat.entropy$entropy)[N])
# abline(v = sort(dat.entropy$entropy)[length(dat.entropy$entropy) - N])
# plot(hc, main = "Clustering: fourier T=24h")

# Plot spectral power
ggplot(dat.mag, aes(x = tissue, y = power)) + geom_bar(stat = "identity") + xlab("Tissue") + ylab("Total 24-h amplitude") + 
  theme_gray(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) + ggtitle("Genome-wide circadian amplitude across tissues")

# Plot heatmaps
PlotRelampHeatmap(M.low, paste0("Low entropy genes.\nN=", length(low.entropy.genes)))
PlotRelampHeatmap(M.high, paste0("High entropy genes.\nN=", length(high.entropy.genes)))

# Plot entropy distributions
plot(density(dat.entropy$entropy[which(!is.na(dat.entropy$entropy))]), 
     main = "Distribution of entropy in rhythmic genes", cex.main = 1.75, cex.lab = 1.75, cex.axis = 1.75)
# draw 2 vertical lines for specifying tissue-spec to tissue-wide
abline(v = sort(dat.entropy$entropy)[N])
abline(v = sort(dat.entropy$entropy)[length(dat.entropy$entropy) - N])

dev.off()

# Figure 1: tissue-specific rhythms ---------------------------------------
# first show MARA results

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)
figcount <- figcount + 1

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long <- LoadActivitiesLong(indir)

act.svd <- GetActSvd(act.long, pval.adj.cutoff = 0.0005)

eigens <- GetEigens(act.svd, comp = 1)
eigens2 <- GetEigens(act.svd, comp = 2)
eigens3 <- GetEigens(act.svd, comp = 3)

# now show some tissue-specific results
# code from find_oscillating_genes.pairs.R
load("Robjs/dat.rhyth.relamp.pvalmin1e-5.pvalmax0.05.relampmax.0.1.meancutoff6.Robj", verbose = T)

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
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# multiplot(eigens$u.plot, eigens$v.plot, n.tissues.rhyth.plot, pairs.plot, layout = jlayout)
multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
multiplot(eigens2$u.plot, eigens2$v.plot, layout = jlayout)
multiplot(eigens3$u.plot, eigens3$v.plot, layout = jlayout)

# jlayout <- matrix(c(1, 2), 2, 1, byrow = TRUE)
# multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
# jlayout <- matrix(c(1, 2), 2, 1, byrow = TRUE)
# multiplot(n.tissues.rhyth.plot, pairs.plot, layout = jlayout)

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
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

hnf4a.exprs <- ggplot(subset(dat.long, gene == "Hnf4a" & experiment == "rnaseq"), aes(x = time, y = exprs)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Hnf4a mRNA expression") + 
  ylab(label = "log2 mRNA expression") +
  xlab("CT") +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

hnf4a.regulated.examp <- ggplot(subset(dat.long, gene == "Upp2" & experiment == "rnaseq"), aes(x = time, y = exprs)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Upp2 (regulated by Bmal1)") + 
  ylab(label = "log2 mRNA expression") +
  xlab("CT") +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))  +
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

mef2c.activity <- ggplot(subset(act.long, gene == "MEF2.A.B.C.D..p2" & experiment == "rnaseq"), aes(x = time, y = exprs)) +
  geom_line() +
  geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
  facet_wrap(~tissue) + 
  xlab("CT") +
  ylab("Activity") + 
  ggtitle("Mef2 motif activity") +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

mef2c.exprs <- ggplot(subset(dat.long, gene == "Mef2c" & experiment == "rnaseq"), aes(x = time, y = exprs)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Mef2c mRNA expression") + 
  ylab(label = "log2 mRNA expression") +
  xlab("CT")  +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

mef2c.regulated.examp <- ggplot(subset(dat.long, gene == "Des" & experiment == "rnaseq"), aes(x = time, y = exprs)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Des (regulated by Mef2c)") + 
  ylab(label = "log2 mRNA expression") +
  xlab("CT")  +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

# mef2c.log2.activity <- ggplot(subset(act.long, gene == "MEF2.A.B.C.D..p2" & experiment == "rnaseq"), aes(x = time, y = 2 ^ exprs - 1)) +
#   geom_line() +
#   # geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
#   facet_wrap(~tissue) + 
#   xlab("CT") +
#   ylab("Activity") + 
#   ggtitle("Mef2 motif activity") +
#   scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
#   theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

mef2c.log2.exprs <- ggplot(subset(dat.long, gene == "Mef2c" & experiment == "rnaseq"), aes(x = time, y = 2 ^ exprs - 1)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Mef2c mRNA expression") + 
  ylab(label = "mRNA expression (linear scale)") +
  xlab("CT")  +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

mef2c.log2.regulated.examp <- ggplot(subset(dat.long, gene == "Des" & experiment == "rnaseq"), aes(x = time, y = 2 ^ exprs - 1)) + 
  geom_line() + 
  facet_wrap(~tissue) +
  ggtitle("Des (regulated by Mef2c)") + 
  ylab(label = "mRNA expression (linear scale)") +
  xlab("CT")  +
  scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
  theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))

# jlayout3 <- matrix(c(1, 2, 3, 4, 5, 6, 0, 7, 8), 3, 3, byrow = TRUE)
jlayout3 <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE)
multiplot(hnf4a.activity, hnf4a.exprs, hnf4a.regulated.examp, 
          mef2c.activity, mef2c.log2.exprs, mef2c.log2.regulated.examp, 
          # mef2c.log2.exprs, mef2c.log2.regulated.examp, 
          layout = jlayout3)

dev.off()

# Figure 4: Alternative promoter ------------------------------------------

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)
figcount <- figcount + 1

load(file = "Robjs/alt_promoter_usage_rel_nr1d1_complex.R", verbose = T)  # tpm.afe, tpm.avg.filt, tpm.fit


# plot examples
hits <- c("Upp2", "Tpm3", "Prkd3", "Insig2", "Ddc")

jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
for (jgene in hits){
  tpm.sub <- subset(tpm.fit, gene_name == jgene)
  jtranscript <- tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1]
  tpm.avg.sub <- subset(tpm.avg.filt, transcript_id == jtranscript)
  exprs.plot <- PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == jgene & experiment == "rnaseq"))
  lin.plot <- ggplot(tpm.avg.sub, aes(y = tpm_norm.avg, x = relamp, label = tissue)) + geom_point() + geom_text() + 
    ggtitle(jgene) + geom_smooth(method = "lm") + ylab("Fractional promoter usage") + xlab("Relative amplitude")
  multiplot(exprs.plot, lin.plot,  layout = jlayout)
}

# plot pval distribution
tpm.summary <- tpm.fit %>%
  group_by(gene_name) %>%
  do(SubsetMinPval(jdf = .))

print(ggplot(tpm.summary, aes(x = pval)) + geom_histogram(binwidth = 0.01) + ggtitle("Distribution of p-values from regression") + xlab("P-value") + ylab("Frequency"))

# load(file = "Robjs/tpm.merged.Robj", verbose = T)
# 
# # for annotating rhythmic genes
# keys <- paste(dat.rhyth.relamp$tissue, dat.rhyth.relamp$gene, sep = ";")
# relampdic <- hash(keys, dat.rhyth.relamp$relamp)
# avgexprsdic <- hash(keys, dat.rhyth.relamp$int.rnaseq)
# 
# tissue.spec <- unique(subset(dat.rhyth.relamp, is.tiss.spec == TRUE)$gene)
# filter.tissues <- c()
# tpm.filt <- subset(tpm.merged, !tissue %in% filter.tissues & gene_name %in% tissue.spec)
# 
# # ~50 seconds
# start.time <- Sys.time()
# tpm.output <- CorrelatePromoterUsageToAmp(tpm.filt, dat.rhyth.relamp, avgexprsdic, filter.tissues, tissue.spec)
# tpm.summary <- tpm.output$tpm.summary
# tpm.avg <- tpm.output$tpm.avg
# tpm.fit <- tpm.output$tpm.fit
# print(Sys.time() - start.time)
# 
# pval.cutoff <- 0.05
# tpm.summary.filt <- subset(tpm.summary, pval <= pval.cutoff)
# sig.hits <- tpm.summary.filt$gene_name
# 
# # plot examples or top hits
# jgene <- "Upp2"
# tpm.sub <- subset(tpm.fit, gene_name == jgene)
# jtranscript <- tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1]
# upp2.altprom.plot <- ggplot(subset(tpm.avg, transcript_id == jtranscript), aes(y = tpm_norm.avg, x = relamp, label = tissue)) + 
#   geom_point() + geom_text() + ggtitle(jgene) + geom_smooth(method = "lm") + xlab("Relative amplitude") + ylab("Fraction promoter usage")
# 
# # plot histogram of p-values
# altprom.histo <- ggplot(tpm.fit, aes(x = pval)) + geom_histogram(binwidth = 0.0075) + xlab("P-value")
# 
# jlayout4 <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# multiplot(altprom.histo, upp2.altprom.plot, layout = jlayout4)

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

exprs.plots <- lapply(jgenes, function(g) return(PlotRnaseqAcrossTissues(subset(dat.long, gene == g & experiment == "rnaseq"), jtitle = paste(g, "RNA-seq profile"))))
biplots <- lapply(jgenes, function(g) return(PlotBiplot(N, g)))

jlayout5 <- matrix(c(1, 4, 2, 5, 3, 6), 3, 2, byrow = TRUE)
do.call(multiplot, c(exprs.plots, biplots, list(layout = jlayout5)))

dev.off()


# Regression model on DHS signal ------------------------------------------

source("scripts/functions/LoadActivitiesLong.R")
pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)

figcount <- figcount + 1

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_DHS/MARA_dhs_gene_region_cutoff_filter/"
zscores.path <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_DHS/MARA_dhs_gene_region_cutoff_filter/"
act.dhs.long <- LoadActivitiesLongDhs(indir, act.file = "Activities", se.file = "StandardError")
zscores <- read.table(zscores.path, header=FALSE, sep="\t", col.names = c("motif", "zscore"))
zscores <- zscores[order(zscores$zscore, decreasing = TRUE), ]

act.sum <- act.dhs.long %>%
  group_by(gene) %>%
  summarise(exprs.sum = sum(exprs))

hits <- c("ONECUT1.2.p2", "NKX2.1.4.p2", "HNF4A_NR2F1.2.p2", "FOXA2.p3", "RFX1..5_RFXANK_RFXAP.p2")
hits.rnaseq <- c("Onecut1", "Nkx2-1", "Hnf4a", "Foxa2", "Rfxank")

# Load RNA-Seq from ENCODE data to plot side by side
load("Robjs/dat.encode.Robj")

dhsplots <- lapply(hits, function(jgene) {
  dat.sub <- SortByTissue(subset(act.dhs.long, gene == jgene), by.var = "exprs")
  return(PlotActivitiesWithSE.dhs(SortByTissue(subset(act.dhs.long, gene == jgene), by.var = "exprs")) + theme_gray(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)))
  })

encodeplots <- lapply(hits.rnaseq, function(jgene) {
  dat.sub <- SortByTissue(subset(dat.encode, gene == jgene & !tissue %in% c("Colon", "Duodenum", "Small Intestine", "Stomach", "Placenta")), by.var = "tpm")
  m <- ggplot(dat.sub, aes(x = tissue, y = tpm)) + geom_bar(stat = "identity") + ggtitle(jgene) + xlab("") + ylab("TPM") +
    theme_gray(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1))
  return(m)
})

for (i in 1:length(hits)){
  dhsgene <- hits[i]
  encodegene <- hits.rnaseq[i]
  dhsplot <- PlotActivitiesWithSE.dhs(subset(act.dhs.long, gene == dhsgene))
  encodeplot <- PlotEncodeRnaseq(subset(dat.encode, gene == encodegene & !tissue %in% c("Colon", "Duodenum", "Small Intestine", "Stomach", "Placenta")))
  
  multiplot(dhsplot, encodeplot, cols = 2)
}

# for (jgene in hits){
#   print(PlotActivitiesWithSE.dhs(subset(act.dhs.long, gene == jgene)) + theme_gray(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)))
# }

dev.off()


# Tars as an example of cooperative binding -------------------------------

pdf(file.path(outdir, paste0("first_year_candidacy_plots_figure", figcount, ".pdf")), height = 8, width = 10.245, paper = "special", onefile = TRUE)
figcount <- figcount + 1

print(PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Tars" & experiment == "rnaseq")) + theme_gray(24))

dev.off()

