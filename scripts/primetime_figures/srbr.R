# 2015-10-20 
# Figures that will be part of paper to analyze this hogenesch dataset

remove.wfat <- TRUE
plot.i <- 1

tissue.order <- c('Liver','BFAT','Kidney','Lung','Adr','Mus','Heart','Aorta','Hypo','Cere','BS')

setwd("/home/yeung/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/srbr"
dir.create(outdir, showWarnings = FALSE)

library(ggplot2)
library(hash)
library(dplyr)
library(reshape2)
library(gplots)
library(penalizedLDA)
library(wordcloud)
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
source(file.path(funcs.dir, "VarianceFunctions.R"))
source(file.path(funcs.dir, "FitRhythmicAcrossPeriods.R"))
source(file.path(funcs.dir, "GetClockGenes.R"))
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/OuterComplex.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/RemoveCommasBraces.R")



# Figure 1C F24 across periods ---------------------------------------------------------------
# Core clock genes
# load(file = "Robjs/dat.fit.scan_periods.Robj")
load(file = "Robjs/dat.fit.periods.genome_wide.min.5_to_30.Robj", verbose=T)
load(file = "Robjs/dat.long.fixed_rik_genes.Robj")

if (remove.wfat){
  dat.long <- subset(dat.long, tissue != "WFAT")
  dat.long$tissue <- factor(dat.long$tissue, levels = unique(dat.long$tissue))
  
  dat.fit.periods.genome_wide.min <- subset(dat.fit.periods.genome_wide.min, tissue != "WFAT")
  dat.fit.periods.genome_wide.min$tissue <- factor(dat.fit.periods.genome_wide.min$tissue, levels = unique(dat.fit.periods.genome_wide.min$tissue))
}
dat.fit.periods.sub <- subset(dat.fit.periods.genome_wide.min, amp > 0.1 & pval < 1e-4)

# order dat.fits by tissue.order
dat.fit.periods.sub$tissue <- factor(dat.fit.periods.sub$tissue, levels = tissue.order)

pdf(file.path(outdir, paste0(plot.i, ".fourier_across_periods.pdf")))
plot.i <- plot.i + 1

xscale_periods <- seq(6, 30, 2)
ggplot(subset(dat.fit.periods.sub), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.sub$period))/55) + geom_vline(xintercept=24, linetype="dotted") + 
  scale_x_continuous(breaks=xscale_periods) + 
  xlab("Period with minimum RSS (h)") + ylab("Number of genes") +
  theme_bw(24) + theme(aspect.ratio=1,
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       axis.text.x=element_text(size=11))
xscale_periods_smaller <- seq(6, 30, 6)
linesize <- 0.1
m2 <- ggplot(subset(dat.fit.periods.sub), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.sub$period))/55) + 
  geom_vline(xintercept=24, linetype="dotted", size = linesize) + 
  geom_vline(xintercept=12, linetype="dotted", size = linesize) +
  scale_x_continuous(breaks=xscale_periods_smaller) + 
  facet_wrap(~tissue) +
  xlab("Period with minimum RSS") + ylab("Count") +
  theme_bw(24) + theme(aspect.ratio=1,
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       axis.text.x=element_text(size=11))
print(m2)

ggplot(dat.fit.periods.sub, aes(x = period, y = ssq.residuals, colour = gene)) + geom_point() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")

dev.off()


# Figure 2 Tissue-modules ----------------------------------------------------------

filt.tiss <- c("WFAT")
# load("Robjs/fits.best.collapsed_models.Robj", verbose=T)
# load("Robjs/fits.best.max_11.collapsed_models.amp_cutoff_0.15.phase_sd.Robj", verbose=T)
# load("Robjs/fits.best.max_11.collapsed_models.amp_cutoff_0.1.phase_sd_diff_avg.Robj", verbose=T)
# load("Robjs/fits.best.max_11.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")
load("Robjs/dat.complex.fixed_rik_genes.Robj")

if (remove.wfat){
  dat.complex <- subset(dat.complex, ! tissue %in% filt.tiss)
}

fits.rhyth <- subset(fits.best, n.params > 0)
fits.rhyth$label <- apply(fits.rhyth, 1, function(row){
  cutoff <- 1
  if (row[8] > cutoff & row[6] > 0){  # amp.avg > cutoff only for n.rhyth > 1
    return(as.character(row[1]))  # return gene
  } else {
    return("")
  }
})

pdf(file.path(outdir, paste0(plot.i, ".tissue_modules.global_stats.pdf")))
plot.i <- plot.i + 1
fits.tspec.sum <- CountModels(subset(fits.best, n.rhyth == 1)) 
m1 <- ggplot(fits.tspec.sum, aes(x = model, y = count)) + geom_bar(stat = "identity") + 
  xlab("") + 
  ylab("Count") + 
  theme_bw(24) + 
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
jmods <- c("Liver", "Adr", "BFAT", "Mus")
m2 <- PlotPolarHistogram(subset(fits.best, n.rhyth == 1 & model %in% jmods), "Count")
multiplot(m1, m2, cols = 2) 

PlotPolarHistogram(subset(fits.best, n.rhyth == 1 & model %in% c(jmods, "Mus")), "Count") + facet_wrap(~tissue, nrow = 1)

dev.off()

pdf(file.path(outdir, paste0(plot.i, ".tissue_modules.tissue_wide.pdf")))
plot.i <- plot.i + 1
dat.mean.rnaseq <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))


# Plot tissue wide genes
fits.tw <- subset(fits.best, n.rhyth >= 8)
genes.tw <- as.character(fits.tw$gene)
#outobj <- PlotHeatmapNconds(fits.tw, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)

s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.tw, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)  

dev.off()

pdf(file.path(outdir, paste0(plot.i, ".overlay_exprs.pdf")))
plot.i <- plot.i + 1

jtiss <- "Adr"
Adr.genes <- as.character(subset(fits.best, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, Adr.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "rnaseq"), Adr.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, Adr.genes, jtiss, avg.method = "median")

jtiss <- "BFAT"
BFAT.genes <- as.character(subset(fits.best, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, BFAT.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "rnaseq"), BFAT.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, BFAT.genes, jtiss, avg.method = "median")

jtiss <- "Mus"
Mus.genes <- as.character(subset(fits.best, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, Mus.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "rnaseq"), Mus.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, Mus.genes, jtiss, avg.method = "median")

jtiss <- "Liver"
Liver.genes <- as.character(subset(fits.best, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, Liver.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "rnaseq"), Liver.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, Liver.genes, jtiss, avg.method = "median")

dev.off()

# Tissue-wide regulation --------------------------------------------------
pdf(file.path(outdir, paste0(plot.i, ".MARA_tissue_wide.pdf")))
plot.i <- plot.i + 1
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long <- LoadActivitiesLong(indir)
act.complex <- TemporalToFrequencyDatLong(act.long, period = 24, n = 8, interval = 6, add.entropy.method = "array")
act.complex$exprs.adj <- act.complex$exprs.transformed * act.complex$frac.weight  # why frac.weight?
act.complex$mod.exprs <- Mod(act.complex$exprs.transformed)
act.complex$mod.exprs.adj <- Mod(act.complex$exprs.adj)

# no WFAT
act.complex <- subset(act.complex, !tissue %in% filt.tiss)

s.act <- SvdOnComplex(act.complex, value.var = "exprs.adj")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
for (comp in seq(1)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 14, pretty.names = TRUE)
  multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
}
dev.off()

# Tissue-specific ----------------------------------------------------

pdf(file.path(outdir, paste0(plot.i, ".MARA_tissue_spec.pdf")))
plot.i <- plot.i + 1
act.long <- subset(act.long, !tissue %in% filt.tiss)
PlotActivitiesWithSE(subset(act.long, gene == "SPIB.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "MEF2.A.B.C.D..p2"))
PlotActivitiesWithSE(subset(act.long, gene == "HNF4A_NR2F1.2.p2"))
dev.off()


# DHS enrichemnt ----------------------------------------------------------

pdf(file.path(outdir, paste0(plot.i, ".dhs_tissue_spec.pdf")))
plot.i <- plot.i + 1
maindir="/home/yeung/projects/tissue-specificity/Robjs/dhs_peaks"
df.out.lst.merged <- lapply(list.files(maindir), function(fname){
  fpath <- file.path(maindir, fname)
  tiss <- strsplit(fname, "\\.")[[1]][[10]]
  print(fpath)
  print(tiss)
  load(fpath)
  df.out.lst.merged$tissue <- tiss
  return(df.out.lst.merged)
})
df.out.lst.merged <- do.call(rbind, df.out.lst.merged)
head(df.out.lst.merged)
jsub <- subset(df.out.lst.merged, gene.type %in% c("jgenes", "jgenes.flat.filt"))
newlab <- hash(as.character(subset(jsub, gene.type == "jgenes")$tissue), paste0(subset(jsub, gene.type == "jgenes")$tissue, "\nN=", subset(jsub, gene.type == "jgenes")$total.genes))
jsub$Tissue <- sapply(as.character(jsub$tissue), function(tiss) newlab[[tiss]])
jsub <- OrderDecreasing(jsub, "Tissue", "frac.n.spec.by.gene")
ggplot(jsub, aes(x = Tissue, y = frac.n.spec.by.gene, fill = gene.type)) + geom_bar(stat = "identity", position = "dodge") + xlab("") + ylab("Fraction of genes with tissue-specific DHS within 10kb") + 
  theme_bw(18) + theme(aspect.ratio = 1, legend.position = "bottom")
dev.off()


# DHS motifs --------------------------------------------------------------

pdf(file.path(outdir, paste0(plot.i, ".penalized_lda_singletons.pdf")))
plot.i <- plot.i + 1
load("Robjs/penalized_lda_liver_matrices.Robj", v=T)

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)

jlambda <- 0.10  # liv only
out.3 <- PenalizedLDA(mat.fgbg.3, labels3, lambda = jlambda, K = 2, standardized = FALSE)
jsize <- sqrt(out.3$discrim[, 1]^2 + out.3$discrim[, 2]^2) * 5 + 0.01
plot(out.3$discrim[, 1], out.3$discrim[, 2], pch = ".", 
     xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others", 
     main="Motifs separating btwn tissues (x-axis) and rhythmicity (y-axis)")
text(out.3$discrim[, 1], out.3$discrim[, 2], names(out.3$x), cex = jsize)
abline(v = 0); abline(h = 0)
# PlotLdaOut2D(out.3, jcex = 0.5)
dev.off()


# WT-KO analysis ----------------------------------------------------------


pdf(file.path(outdir, paste0(plot.i, ".wtko.pdf")))
plot.i <- plot.i + 1
load("Robjs/wtko.dat.wtko.hog.commongenes.Robj", v=T)
load("Robjs/wtko.fits.all.long.commongenes.Robj", v=T)
# load("Robjs/wtko.dat.complex.wtko.Robj", v=T)
load("Robjs/wtko.dat.complex.wtko.hog.Robj", v=T)

genes.liv <- as.character(subset(fits.best, model == "Liver")$gene)
genes.tw <- as.character(subset(fits.best, n.rhyth >= 8)$gene)

fits.sub.liv <- subset(fits.all.long.wtkohog, gene %in% genes.liv & model != "")

# genes.liv.wtliv <- as.character(subset(fits.sub.liv, model %in% c("WT,Liver", "WT;KO;Liver", "WT;Liver", "KO;WT,Liver", "Liver", "WT", "WT,Liver", "WT;Liver"))$gene)
# genes.liv.wtliv <- as.character(subset(fits.sub.liv, model %in% c("WT,Liver", "WT;KO;Liver", "WT;Liver", "KO;WT,Liver"))$gene)
# genes.liv.wtliv <- as.character(subset(fits.sub.liv, model %in% c("WT,Liver", "WT;Liver", "Liver", "WT", "WT,Liver", "WT;Liver"))$gene)
genes.liv.wtliv <- as.character(subset(fits.sub.liv, model %in% c("WT,Liver", "WT;Liver", "Liver", "WT"))$gene)
print(paste("N genes:", length(genes.liv.wtliv)))

# s.liv.wtliv <- SvdOnComplex(subset(dat.complex.wtko, gene %in% genes.liv.wtliv), value.var = "exprs.transformed")
s.liv.wtliv <- SvdOnComplex(subset(dat.complex.wtko.hog, gene %in% genes.liv.wtliv), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.liv.wtliv, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)
dev.off()


# WT-KO genes motif enrichment --------------------------------------------


pdf(file.path(outdir, paste0(plot.i, ".wtko.pdf")))
plot.i <- plot.i + 1

mat.fgbg.lab.lst.3 <- SetUpMatForLda(subset(mat.fg, gene %in% genes.liv.wtliv), mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)


jlambda <- 0.1  # liv only
out.3 <- PenalizedLDA(mat.fgbg.3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

jsize <- sqrt(out.3$discrim[, 1]^2 + out.3$discrim[, 2]^2) * 5 + 0.01
plot(out.3$discrim[, 1], out.3$discrim[, 2], pch = ".", 
     xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others", 
     main="Motifs separating btwn tissues (x-axis) and rhythmicity (y-axis)")
text(out.3$discrim[, 1], out.3$discrim[, 2], names(out.3$x), cex = jsize)
abline(v = 0); abline(h = 0)

# PlotLdaOut2D(out.3, jcex = 0.5)

dev.off()
