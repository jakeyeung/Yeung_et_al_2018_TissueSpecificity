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
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/FourierFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/LongToMat.R")



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
  xlab("Period with best fit (h)") + ylab("Number of genes") +
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
  xlab("Period with best fit (h)") + ylab("Count") +
  theme_bw(24) + theme(aspect.ratio=1,
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       axis.text.x=element_text(size=11))
print(m2)

ggplot(dat.fit.periods.sub, aes(x = period, y = ssq.residuals, colour = gene)) + geom_point() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")

dev.off()


# Genome-wide amplitudes --------------------------------------------------

# from fourier directory total_variance.noise_floor.R
pdf(file.path(outdir, paste0(plot.i, ".genomewide_amplitude.pdf")))
plot.i <- plot.i + 1

load("Robjs/dat.fit.Robj", v=T); dat.fit.24 <- dat.fit
dat.fit.24 <- subset(dat.fit.24, tissue != "WFAT")
dat.fit.24 <- dat.fit.24[order(dat.fit.24$amp, decreasing = TRUE), ]

amp.thres <- seq(from = 0, to = max(dat.fit.24$amp), by = 0.15)

pval.cutoff <- 0.01
dat.fit.24.ngenes.thres <- subset(dat.fit.24, pval < pval.cutoff) %>%
  group_by(tissue) %>%
  do(NGenesByAmp.long(., amp.thres))
dat.fit.24.ngenes.thres$rhyth <- as.factor(24)

# order by total genes
ngenes.sum <- dat.fit.24.ngenes.thres %>%
  group_by(tissue) %>%
  summarise(total = sum(n.genes)) %>%
  arrange(desc(total))
dat.fit.24.ngenes.thres$tissue <- factor(as.character(dat.fit.24.ngenes.thres$tissue), levels = ngenes.sum$tissue)
ggplot(subset(dat.fit.24.ngenes.thres, rhyth == 24), aes(x = 2 * amp.thres, y = n.genes, colour = tissue)) + geom_line() + 
  theme_bw(24) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank()) +
  xlab("Log2 Fold Change") + ylab("# Genes") + xlim(c(0, 5)) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000)) + 
  geom_vline(xintercept = 1, linetype = "dotted")

dev.off()


# Examples to motivate ----------------------------------------------------

pdf(file.path(outdir, paste0(plot.i, ".gene_exprs_examps.pdf")))
plot.i <- plot.i + 1

# load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj")

library(hash)

jgenes <- c("Nr1d1", "Dbp", "Ndrg1", "Pi4k2a", "Slc44a1")
jgenes <- rev(jgenes)
m.list <- list()
i <- 1
# jexperiments <- c("array", "rnaseq")
# for (jexp in jexperiments){
jexp <- "array"
  for (jgene in jgenes){
    m <- PlotGeneByRhythmicParameters(fits.best, subset(dat.long, experiment == jexp), 
                                      jgene, amp.filt = 0.2, jtitle=jgene, facet.rows = 1, jcex = 8,
                                      pointsize = 0)
    # m <- PlotGeneByRhythmicParameters(fits.best, subset(dat.long), jgene, amp.filt = 0.2, jtitle=jgene, facet.rows = 1, jcex = 8)
    print(m)
    m.list[[i]] <- m
    i <- i + 1
  }
do.call(multiplot, m.list)
# }
# multiplot(m.list[[1]], m.list[[2]], m.list[[3]], m.list[[4]], cols = 2)
dev.off()


pdf(file.path(outdir, paste0(plot.i, ".model_selection_summary.pdf")))
plot.i <- plot.i + 1
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")

fits.best$n.rhyth.fac <- as.factor(sapply(as.numeric(fits.best$n.rhyth), function(n) NrhythToStr(n)))
ggplot(subset(fits.best, n.rhyth != 0), aes(x = as.factor(n.rhyth.fac), y = 2 * amp.avg)) + geom_boxplot() + xlab("# rhythmic tissues") + ylab("Log fold change")  + theme_bw(16) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

# Figure 2 Tissue-modules ----------------------------------------------------------

pdf(file.path(outdir, paste0(plot.i, ".tissue_modules.global_stats.pdf")))
plot.i <- plot.i + 1

filt.tiss <- c("WFAT")
# load("Robjs/fits.best.collapsed_models.Robj", verbose=T)
# load("Robjs/fits.best.max_11.collapsed_models.amp_cutoff_0.15.phase_sd.Robj", verbose=T)
# load("Robjs/fits.best.max_11.collapsed_models.amp_cutoff_0.1.phase_sd_diff_avg.Robj", verbose=T)
# load("Robjs/fits.best.max_11.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")
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

# count based on amp
amp.thres <- seq(from = 0, to = max(dat.fit.24$amp), by = 0.15)

fits.best$n.rhyth.lab <- sapply(fits.best$n.rhyth, function(n){
  if (n >= 8){
    return("8-11")
  } else if (n == 1){
    return("1")
  } else if (n <= 3){
    return("2-4")
  } else {
    return("5-7")
  }
})
fits.counts.by.amp <- subset(fits.best, n.rhyth > 0) %>%
  group_by(n.rhyth.lab) %>%
  do(NGenesByAmp.long(., amp.thres, labelid = "n.rhyth.lab", varid = "amp.avg", outlabel = "n.rhyth.lab"))
ggplot(fits.counts.by.amp, aes(x = 2 * amp.thres, y = n.genes, group = n.rhyth.lab, colour = as.factor(n.rhyth.lab))) + geom_line() + 
  geom_line(size = 2) + 
  theme_bw(20) +
  labs(colour = "# Rhythmic\nTissues") + 
  theme(aspect.ratio=1, 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Avg Amplitude of Rhythmic Tissues") + ylab("# Genes") + xlim(c(0.15, 6)) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000)) + 
  geom_vline(xintercept = 2.8, linetype = "dotted") + 
  scale_colour_brewer(palette = "Spectral")

fits.tspec.sum <- CountModels(subset(fits.best, n.rhyth == 1)) 
m1 <- ggplot(fits.tspec.sum, aes(x = model, y = count)) + geom_bar(stat = "identity") + 
  xlab("") + 
  ylab("Count") + 
  theme_bw(24) + 
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
jmods <- c("Liver", "Adr", "BFAT", "Mus")
jmods2 <- c("Adr", "BFAT", "Mus")
m2 <- PlotPolarHistogram(subset(fits.best, n.rhyth == 1 & model %in% jmods), "Count")
m3 <- PlotPolarHistogram(subset(fits.best, n.rhyth == 1 & model %in% jmods2), "Count")
multiplot(m1, m2, cols = 2) 
multiplot(m1, m3, cols = 2) 

PlotPolarHistogram(subset(fits.best, n.rhyth == 1 & model %in% c(jmods, "Mus")), "Count") + facet_wrap(~tissue, nrow = 1)
PlotPolarHistogram(subset(fits.best, n.rhyth == 1 & model %in% c(jmods2)), "Count") + facet_wrap(~tissue, nrow = 1)

dev.off()

dat.mean.rnaseq <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))


pdf(file.path(outdir, paste0(plot.i, ".tissue_modules.tissue_wide.pdf")))
plot.i <- plot.i + 1

# Plot tissue wide genes
fits.tw <- subset(fits.best, n.rhyth >= 8)
genes.tw <- as.character(fits.tw$gene)
#outobj <- PlotHeatmapNconds(fits.tw, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)

s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.tw, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough=TRUE, ylab="Phase (CT)", jtitle = "")
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.tw$v.plot, eigens.tw$u.plot, layout = jlayout)  

dev.off()

pdf(file.path(outdir, paste0(plot.i, ".overlay_exprs.pdf")), paper = "USr")
plot.i <- plot.i + 1

# plot all 4 tissues
jtissues <- c("Liver", "Adr", "BFAT", "Mus")
jtissues2 <- c("Adr", "BFAT", "Mus")
# jtissues <- c("Adr")
genes.lst <- list()
for (tiss in jtissues){
  genes.lst[[tiss]] <- as.character(subset(fits.best, model == tiss)$gene)
}
jgenes <- unlist(genes.lst)
dat.sub <- subset(dat.long, experiment == "array" & gene %in% jgenes & tissue %in% jtissues)
# rearrange tissues
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = jtissues)
m <- PlotOverlayTimeSeries(dat.sub, genes.lst, tissues = jtissues, jalpha = 0.025, jtitle = paste0(""))
m2 <- PlotOverlayTimeSeries(dat.sub, genes.lst, tissues = jtissues2, jalpha = 0.025, jtitle = paste0(""))
print(m)
print(m2)

# fits.sub <- subset(fits.best, model %in% jtissues)
# models.hash <- hash(as.character(fits.sub$gene), as.character(fits.sub$model))
# dat.mean.rnaseq.sub <- subset(dat.mean.rnaseq, gene %in% as.character(fits.sub$gene))
# dat.mean.rnaseq.sub$model <- sapply(as.character(dat.mean.rnaseq.sub$gene), function(g) models.hash[[g]])
# m2 <- PlotMeanExprsOfModel(dat.mean.rnaseq.sub, unlist(genes.lst, use.names = FALSE), jtissues) + facet_wrap(~model, nrow = 1)
# print(m2)

plots.lst <- list()
for (tiss in jtissues){
  genes <- genes.lst[[tiss]]
  plots.lst[[tiss]] <- PlotMeanExprsOfModel(dat.mean.rnaseq, genes, tiss, avg.method = "median") + theme_bw(6)
}
multiplot(plots.lst[[1]], plots.lst[[2]], plots.lst[[3]], plots.lst[[4]], cols = 4)

jtiss <- "Adr"
Adr.genes <- as.character(subset(fits.best, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, Adr.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "array"), Adr.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, Adr.genes, jtiss, avg.method = "median")

jtiss <- "BFAT"
BFAT.genes <- as.character(subset(fits.best, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, BFAT.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "array"), BFAT.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, BFAT.genes, jtiss, avg.method = "mean")

# plot mean exprs of BFAT with Mef2c accumulation
jsub.mean <- subset(dat.mean.rnaseq, gene %in% BFAT.genes)
# annotate mef2c mrna accumulation
jsub.mef2c <- subset(dat.mean.rnaseq, gene == "Mef2c")
jsub.hash <- hash(as.character(jsub.mef2c$tissue), jsub.mef2c$exprs.mean)
jsub.mean$mef2c <- sapply(as.character(jsub.mean$tissue), function(tiss) jsub.hash[[tiss]])

jsub.mean.meds <- jsub.mean %>%
  group_by(tissue) %>%
  summarise(exprs.med = median(exprs.mean),
            mef2c = unique(mef2c),
            upper = quantile(exprs.mean, probs = 0.75),
            lower = quantile(exprs.mean, probs = 0.25))
ggplot(jsub.mean.meds, aes(x = mef2c, y = exprs.med, label = tissue)) + 
  geom_point(size = 2) + 
  geom_text(aes(x = mef2c, y = upper + 0.3), size = 7.5) + 
  xlab("Mean Mef2c mRNA accumulation") + 
  ylab("Mean exprs of gene") + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.25)) + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# plot in a way that clusters other genes
jgene <- "Mef2c"
m1 <- PlotGeneByRhythmicParameters(fits.best, subset(dat.long, experiment == "array"), 
                                   jgene, amp.filt = 0.2, jtitle=jgene, facet.rows = 1, jcex = 8,
                                   pointsize = 0) + theme_bw(24) + theme(legend.position = "none", aspect.ratio = 1)
print(m1)


jtiss <- "Mus"
Mus.genes <- as.character(subset(fits.best, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, Mus.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "array"), Mus.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, Mus.genes, jtiss, avg.method = "median")

jtiss <- "Liver"
Liver.genes <- as.character(subset(fits.best, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, Liver.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "array"), Liver.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, Liver.genes, jtiss, avg.method = "median")

dev.off()

# Tissue-wide regulation --------------------------------------------------
pdf(file.path(outdir, paste0(plot.i, ".MARA_tissue_wide.pdf")))
plot.i <- plot.i + 1
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"

act.long <- LoadActivitiesLong(indir)
act.complex <- TemporalToFrequencyDatLong(act.long, period = 24, n = 8, interval = 6, add.entropy.method = "array")

# no WFAT
act.complex <- subset(act.complex, !tissue %in% filt.tiss)

s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
for (comp in seq(1)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 25, pretty.names = TRUE, peak.to.trough = TRUE, jtitle = "")
  multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
}

# tissue wide genes
twmaradir <- "/home/yeung/projects/tissue-specificity/results/MARA/bic_modules/TissueWide.centeredTRUE"

act.long <- LoadActivitiesLong(twmaradir)
act.complex <- TemporalToFrequencyDatLong(act.long, period = 24, n = 8, interval = 6, add.entropy.method = "array")

# no WFAT
act.complex <- subset(act.complex, !tissue %in% filt.tiss)

s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
for (comp in seq(1)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 25, pretty.names = TRUE, peak.to.trough = TRUE, jtitle = "")
  multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
}

# Expression of HIC1 across tissues
print(PlotGeneAcrossTissues(subset(dat.long, gene == "Hic1" & experiment == "rnaseq"), make.pretty = TRUE))


dev.off()

# Tissue-specific ----------------------------------------------------

pdf(file.path(outdir, paste0(plot.i, ".MARA_tissue_spec.pdf")))
plot.i <- plot.i + 1
act.long <- subset(act.long, !tissue %in% filt.tiss)
jmotifs.lst <- list("BFAT"="MEF2.A.B.C.D..p2", "Adr"="HNF4A_NR2F1.2.p2")
for (jtiss in names(jmotifs.lst)){
  jmotif <- jmotifs.lst[[jtiss]]
  m <- PlotActivitiesWithSE(subset(act.long, gene == jmotif & tissue == jtiss & experiment == "array")) + theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
  m1 <- PlotActivitiesWithSE(subset(act.long, gene == jmotif & experiment == "array")) + theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
  print(m)
  print(m1)
}

jmotif <- "MEF2.A.B.C.D..p2"
act.sub <- subset(act.long, gene == jmotif & experiment == "array")
act.sub.mean <- act.sub %>%
  group_by(tissue) %>%
  summarise(mean.exprs = mean(exprs)) %>%
  arrange(desc(mean.exprs))
# act.sub$tissue <- factor(as.character(act.sub$tissue), levels = c("Mus", "Heart", "BS", "Hypo", "Cere", "Aorta", "Lung", "Adr", "BFAT", "Kidney", "Liver"))
act.sub$tissue <- factor(as.character(act.sub$tissue), levels = as.character(act.sub.mean$tissue))
PlotActivitiesWithSE(act.sub) + theme_bw() + theme(aspect.ratio = 1, legend.position = "none")

dat.sub <- subset(dat.long, gene == "Mef2c" & experiment == "array")
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = as.character(act.sub.mean$tissue))
PlotGeneAcrossTissues(dat.sub) + theme_bw()  + theme(aspect.ratio = 1, legend.position = "none")

# plot mean exprs with mRNA accumulation of mef2c


PlotActivitiesWithSE(subset(act.long, gene == "SPIB.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "HNF4A_NR2F1.2.p2"))
dev.off()


# Mean expression of modules ----------------------------------------------

# PlotMeanExprsOfModel(dat.mean, genes, jmodel, sorted = TRUE, avg.method = "mean")

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
jsub <- subset(df.out.lst.merged, gene.type %in% c("jgenes"))
newlab <- hash(as.character(subset(jsub, gene.type == "jgenes")$tissue), paste0(subset(jsub, gene.type == "jgenes")$tissue, "\nN=", subset(jsub, gene.type == "jgenes")$total.genes))
jsub$Tissue <- sapply(as.character(jsub$tissue), function(tiss) newlab[[tiss]])
jsub <- OrderDecreasing(jsub, "Tissue", "frac.n.spec.by.gene")
ggplot(jsub, aes(x = Tissue, y = frac.n.spec.by.gene, fill = gene.type)) + geom_bar(stat = "identity", position = "dodge") + xlab("") + ylab("Fraction of genes with tissue-specific DHS within 10kb") + 
  theme_bw(18) + theme(aspect.ratio = 1, legend.position = "none")
dev.off()


# DHS motifs --------------------------------------------------------------

pdf(file.path(outdir, paste0(plot.i, ".penalized_lda_singletons.pdf")))
plot.i <- plot.i + 1
# load("Robjs/penalized_lda_liver_matrices.Robj", v=T)
load("Robjs/penalized_lda_matrices.dist10kb.jcutoff2.Robj", v=T)

# WITH LIVER vs FLAT
mat.fgbg.lab.lst <- SetUpMatForLda(mat.fg, mat.bg, has.peaks = TRUE)
mat.fgbg <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)
jlambda <- 0.1  # liv only
out <- PenalizedLDA(mat.fgbg, labels, lambda = jlambda, K = 1, standardized = FALSE)
jsize <- sqrt(out$discrim[, 1]^2 * 5) + 0.01
PlotLdaOut(out, jcex = 0.6, jtitle = "")


# WITH LIVER vs FLAT vs NONLIVER
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

# show an example
jgene <- "Mreg"
jsub <- subset(dat.wtko.hog, gene == jgene & experiment != "array")
jsub <- jsub %>%
  group_by(tissue) %>%
  mutate(exprs.center = exprs - mean(exprs))
jsub$exprs <- jsub$exprs.center + 10

jsub$tissue <- factor(jsub$tissue, levels = c("Liver", "WT", "KO"))
PlotGeneAcrossTissues(jsub) + theme_bw() + theme(aspect.ratio = 1)

# s.liv.wtliv <- SvdOnComplex(subset(dat.complex.wtko, gene %in% genes.liv.wtliv), value.var = "exprs.transformed")
s.liv.wtliv <- SvdOnComplex(subset(dat.complex.wtko.hog, gene %in% genes.liv.wtliv), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.liv.wtliv, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.tw$v.plot, eigens.tw$u.plot + xlim(0, 2), layout = jlayout)
multiplot(eigens.tw$v.plot, eigens.tw$u.plot, layout = jlayout)

# plot heatmap
# library(gplots)
# fits.best.sub <- subset(fits.best, gene %in% genes.liv.wtliv)
# dat.sub <- subset(dat.wtko.hog, gene %in% genes.liv.wtliv)
# dat.sub$experiment <- "array"
# PlotHeatmapNconds(fits.best.sub, dat.sub, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)

dev.off()


# WT-KO genes motif enrichment --------------------------------------------


pdf(file.path(outdir, paste0(plot.i, ".motifs.pdf")))
plot.i <- plot.i + 1

mat.fgbg.lab.lst <- SetUpMatForLda(subset(mat.fg, gene %in% genes.liv), mat.bg, has.peaks = TRUE)
mat.fgbg <- mat.fgbg.lab.lst$mat.fgbg; labels <- mat.fgbg.lab.lst$labels
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)
jlambda <- 0.1  # liv only
out <- PenalizedLDA(mat.fgbg, labels, lambda = jlambda, K = 1, standardized = FALSE)
PlotLdaOut(out, jcex = 0.5)

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