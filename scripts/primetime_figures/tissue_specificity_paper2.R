# 2016-06-02
# Clean up figures


# Inits -------------------------------------------------------------------

remove.wfat <- TRUE
plot.i <- 1

tissue.order <- c('Liver','BFAT','Kidney','Lung','Adr','Mus','Heart','Aorta','Hypo','Cere','BS')

setwd("/home/yeung/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/primetime_paper_redo"
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



# Fourier Analysis --------------------------------------------------------

# Fourier across tissues
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


xscale_periods <- seq(6, 30, 2)
plot.periods.all <- ggplot(subset(dat.fit.periods.sub), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.sub$period))/55) + geom_vline(xintercept=24, linetype="dotted") + 
  scale_x_continuous(breaks=xscale_periods) + 
  xlab("Period with best fit [h]") + ylab("Number of genes") +
  theme_bw(24) + theme(aspect.ratio=1,
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       axis.text.x=element_text(size=11))
xscale_periods_smaller <- seq(6, 30, 6)
linesize <- 0.1
plot.periods.bytiss <- ggplot(subset(dat.fit.periods.sub), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.sub$period))/55) + 
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

### BEGIN FOURIER ANALYSIS ###
load("Robjs/dat.complex.all_T.rbinded.Robj", verbose=T)

if (remove.wfat){
  dat.complex.all_T <- subset(dat.complex.all_T, tissue != "WFAT")
}

dat.var.s <- dat.complex.all_T %>%
  group_by(period, tissue) %>%
  summarise(sum_sqr_mod = sum(Mod(exprs.transformed) ^ 2) * 2) %>%  # * 2 to consider symmetrical frequencies
  mutate(period.factor = signif(period, digits = 3))
dat.var.s$period.factor <- factor(dat.var.s$period.factor, 
                                  levels = sort(unique(dat.var.s$period.factor), decreasing = TRUE))
dat.var.all <- dat.var.s %>%
  group_by(tissue) %>%
  summarise(sum_sqr_mod.total = sum(sum_sqr_mod) * 2)  # * 2 considers symmetrical frequencies
dat.var.all$tissue <- factor(dat.var.all$tissue,
                             levels = dat.var.all$tissue[order(dat.var.all$sum_sqr_mod.total, decreasing = TRUE)])

# period.factor: everything not 24 and 12 should be called Other
dat.var.s$period.factor.cond <- sapply(dat.var.s$period.factor, function(f){
  f <- as.character(f)
  if (f == "24" | f == "12"){
    return(as.factor(f))
  }
  else{
    return(as.factor("Other"))
  }
})
dat.var.s <- dat.var.s %>%
  group_by(tissue) %>%
  arrange(period.factor.cond)

dat.var.s1_adj <- dat.var.s %>%
  group_by(tissue) %>%
  mutate(s1_normalized = sum_sqr_mod / sum(sum_sqr_mod))

dat.var.s1_adj.24 <- subset(dat.var.s1_adj, period.factor.cond == "24")
dat.var.s1_adj.24$tissue <- factor(dat.var.s1_adj.24$tissue,
                                   levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])

dat.var.s1_adj.12 <- subset(dat.var.s1_adj, period.factor.cond == "12")
dat.var.s1_adj.12$tissue <- factor(dat.var.s1_adj.12$tissue,
                                   levels = dat.var.s1_adj.12$tissue[order(dat.var.s1_adj.12$s1_normalized, decreasing = TRUE)])

# for ordering the facet_wrap plot across tissues by total variance
dat.var.s$tissue <- factor(dat.var.s$tissue,
                           levels = dat.var.all$tissue[order(dat.var.all$sum_sqr_mod.total, decreasing = TRUE)])

# for ordering the facet_wrap plot across tissues by 24h variance
dat.var.s$tissue <- factor(dat.var.s$tissue,
                           levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$sum_sqr_mod, decreasing = TRUE)])
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
barwidth=0.7
plot.spectral.power <- ggplot() + 
  geom_bar(data=subset(dat.var.s, period.factor.cond != "Other"), mapping=aes(x = tissue, y = sum_sqr_mod, fill = period.factor.cond),
           stat="identity", position=position_dodge(), colour="black", width = barwidth) +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       aspect.ratio = 0.5, 
                       legend.position="bottom") + xlab("") + ylab("Spectral power [log2]^2") +
  scale_fill_manual(name = "Period [h]", drop=TRUE, values=cbPalette)
# add noise floor
periods <- sort(unique(dat.var.s$period))
noise.components <- periods[which(24 %% periods != 0)]

# mean for all tissues
dat.var.s.noisefloor <- subset(dat.var.s, period %in% noise.components) %>%
  group_by(period.factor) %>%
  summarise(sum_sqr_mod.mean = mean(sum_sqr_mod), 
            sum_sqr_mod.var = var(sum_sqr_mod), 
            sum_sqr_mod.min = min(sum_sqr_mod), 
            sum_sqr_mod.max = max(sum_sqr_mod))
# tissue by tissue
dat.var.s.noisefloor.bytiss <- subset(dat.var.s, period %in% noise.components) %>%
  group_by(tissue) %>%
  summarise(sum_sqr_mod.mean = mean(sum_sqr_mod))
lstart <- 0.65
dat.var.s.noisefloor.bytiss$start <- seq(lstart, by = 1, length.out = nrow(dat.var.s.noisefloor.bytiss))
dat.var.s.noisefloor.bytiss$end <- seq(lstart + barwidth, by = 1, length.out = nrow(dat.var.s.noisefloor.bytiss))

# add noise floor that is tissuewide
noise.floor <- mean(dat.var.s.noisefloor$sum_sqr_mod.mean)
plot.spectral.power.noiseflr <- plot.spectral.power + geom_hline(aes(yintercept = noise.floor), linetype="dotted")

# add tissue-specific noise floor
plot.spectral.power.noiseflr.bytiss <- plot.spectral.power + geom_segment(data = dat.var.s.noisefloor.bytiss, mapping=aes(x=start, xend=end, y=sum_sqr_mod.mean, yend=sum_sqr_mod.mean), linetype = "dotted")


# plot normalized spectral power
dat.var.s1_adj$tissue <- factor(dat.var.s1_adj$tissue,
                                levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])
# add tissue-specific noise floor
dat.var.s1_adj.noisefloor.bytiss <- subset(dat.var.s1_adj, period %in% noise.components) %>%
  group_by(tissue) %>%
  summarise(s1_normalized = mean(s1_normalized))
dat.var.s1_adj.noisefloor.bytiss$start <- seq(lstart, by = 1, length.out = nrow(dat.var.s1_adj.noisefloor.bytiss))
dat.var.s1_adj.noisefloor.bytiss$end <- seq(lstart + barwidth, by = 1, length.out = nrow(dat.var.s1_adj.noisefloor.bytiss))

plot.normalized.spectral.power <- ggplot() + 
  geom_bar(data=subset(dat.var.s1_adj, period.factor.cond != "Other"), mapping=aes(x = tissue, y = s1_normalized, fill = period.factor.cond),
           stat="identity", position=position_dodge(), colour="black", width = barwidth) +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       aspect.ratio = 0.5, 
                       legend.position="bottom") + xlab("") + ylab("Normalized spectral power") +
  scale_fill_manual(name = "Period [h]", drop=TRUE, values=cbPalette)
plot.normalized.spectral.power <- plot.normalized.spectral.power + geom_segment(data = dat.var.s1_adj.noisefloor.bytiss, mapping=aes(x=start, xend=end, y=s1_normalized, yend = s1_normalized), linetype = "dotted")
### END FOURIER ANALYSIS ###

### BEGIN GENOMEWIDE AMP ANALYSIS ###
load("Robjs/dat.fit.Robj", v=T); dat.fit.24 <- dat.fit
load("Robjs/dat.fit.period.12.Robj", v=T); dat.fit.12 <- dat.fit
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)

source("scripts/functions/PlotGeneAcrossTissues.R", v=T)

if (remove.wfat){
  dat.fit.24 <- subset(dat.fit.24, tissue != "WFAT")
  dat.fit.12 <- subset(dat.fit.12, tissue != "WFAT")
}

dat.fit.12 <- dat.fit.12[order(dat.fit.12$amp, decreasing = TRUE), ]
dat.fit.24 <- dat.fit.24[order(dat.fit.24$amp, decreasing = TRUE), ]

amp.thres <- seq(from = 0, to = max(dat.fit.12$amp, dat.fit.24$amp), by = 0.15)

pval.cutoff <- 0.01
dat.fit.24.ngenes.thres <- subset(dat.fit.24, pval < pval.cutoff) %>%
  group_by(tissue) %>%
  do(NGenesByAmp.long(., amp.thres))
dat.fit.12.ngenes.thres <- subset(dat.fit.12, pval < pval.cutoff) %>%
  group_by(tissue) %>%
  do(NGenesByAmp.long(., amp.thres))

# combine the two
dat.fit.24.ngenes.thres$rhyth <- as.factor(24)
dat.fit.12.ngenes.thres$rhyth <- as.factor(12)
dat.fit.ngenes.thres <- rbind(dat.fit.24.ngenes.thres, dat.fit.12.ngenes.thres)

# plot 24h and 12h rhythms for each tissue
# ggplot(dat.fit.ngenes.thres, aes(x = 2 * amp.thres, y = n.genes + 1, colour = rhyth)) + geom_line() + 
#   facet_wrap(~tissue) + theme_bw(24) + ggtitle(paste("12 and 24 hour rhythms. Pval >", pval.cutoff)) +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   xlab("Amplitude (fold change peak/trough)") + ylab("Number of Genes") + xlim(c(0, 5)) + 
#   scale_y_log10()

# order by total genes
ngenes.sum <- dat.fit.ngenes.thres %>%
  group_by(tissue) %>%
  summarise(total = sum(n.genes)) %>%
  arrange(desc(total))
dat.fit.ngenes.thres$tissue <- factor(as.character(dat.fit.ngenes.thres$tissue), levels = ngenes.sum$tissue)
plot.genomewide.amps <- ggplot(subset(dat.fit.ngenes.thres, rhyth == 24), aes(x = 2 * amp.thres, y = n.genes, colour = tissue)) + geom_line(size = 2.5) + 
  theme_bw(24) +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  xlab("Log2 Fold Change") + ylab("# Genes") + xlim(c(0, 7)) + 
  scale_y_log10(breaks = c(10, 100, 1000)) + 
  geom_vline(xintercept = 2, linetype = "dotted") + 
  scale_colour_brewer(palette = "Spectral")
### END GENOMEWIDE AMP ANALYSIS ###


# print output
pdf(file.path(outdir, paste0(plot.i, ".overview_of_dataset.pdf")))
plot.i <- plot.i + 1
print(plot.periods.all)
print(plot.periods.bytiss)
print(plot.spectral.power.noiseflr)
print(plot.spectral.power.noiseflr.bytiss)
print(plot.normalized.spectral.power)
print(plot.genomewide.amps)
dev.off()


# Nconds Summary ----------------------------------------------------------

### BEGIN AMPS OF MODULES ###
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

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
plots.nconds.amps <- ggplot(fits.counts.by.amp, aes(x = 2 * amp.thres, y = n.genes, group = n.rhyth.lab, colour = as.factor(n.rhyth.lab))) + geom_line() + 
  geom_line(size = 2) + 
  theme_bw(20) +
  labs(colour = "# Rhythmic\nTissues") + 
  theme(aspect.ratio=0.5, 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Avg Amplitude of Rhythmic Tissues") + ylab("# Genes") + xlim(c(0.15, 6)) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000)) + 
  geom_vline(xintercept = 2.8, linetype = "dotted") + 
  scale_colour_brewer(palette = "Spectral")

### END AMPS OF MODULES ###

### BEGIN EXAMPLES ###

jgenes <- c("Dbp", "Ndrg1", "Pi4k2a", "Slc44a1")
jgenes <- rev(jgenes)
plots.nconds.examples <- list()
i <- 1
# jexperiments <- c("array", "rnaseq")
# for (jexp in jexperiments){
jexp <- "array"
for (jgene in jgenes){
  m <- PlotGeneByRhythmicParameters(fits.best, subset(dat.long, experiment == jexp & tissue != "WFAT"), 
                                    jgene, amp.filt = 0.2, jtitle=jgene, facet.rows = 1, jcex = 8,
                                    pointsize = 0)
  plots.nconds.examples[[i]] <- m
  i <- i + 1
}
### END EXAMPLES ###


### BEGIN SUMMARY OF MODULES ###

# Tissue-specific
fits.tspec.sum <- CountModels(subset(fits.best, n.rhyth == 1)) 
plots.ts.counts <- ggplot(fits.tspec.sum, aes(x = model, y = count)) + geom_bar(stat = "identity") + 
  xlab("") + 
  ylab("Count") + 
  theme_bw(20) + 
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
jmods <- c("Liver", "Adr", "BFAT", "Mus")
plots.phase.histo <- PlotPolarHistogram(subset(fits.best, n.rhyth == 1 & model %in% jmods), "Count") + theme(strip.text.x = element_text(size = 16))

# tissue-wide
load("Robjs/dat.complex.fixed_rik_genes.Robj")
if (remove.wfat){
  dat.complex <- subset(dat.complex, tissue != "WFAT")
}
fits.tw <- subset(fits.best, n.rhyth >= 8)
genes.tw <- as.character(fits.tw$gene)

s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.tw, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 5, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

### END SUMMARY OF MODULES

# print output
pdf(file.path(outdir, paste0(plot.i, ".nconds_summary.pdf")))
plot.i <- plot.i + 1

print(plots.nconds.amps)
do.call(multiplot, plots.nconds.examples)
multiplot(plots.ts.counts, plots.phase.histo)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)  

dev.off()

# Regulation by MARA analysis ---------------------------------------------

### BEGIN TISSUEWIDE REGULATORS ###
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long <- LoadActivitiesLong(indir)
act.complex <- TemporalToFrequencyDatLong(act.long, period = 24, n = 8, interval = 6, add.entropy.method = "array")
act.complex$exprs.adj <- act.complex$exprs.transformed * act.complex$frac.weight  # why frac.weight?
act.complex$mod.exprs <- Mod(act.complex$exprs.transformed)
act.complex$mod.exprs.adj <- Mod(act.complex$exprs.adj)

# no WFAT
act.complex <- subset(act.complex, tissue != "WFAT")

s.act <- SvdOnComplex(act.complex, value.var = "exprs.adj")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
for (comp in seq(1)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 14, pretty.names = TRUE)
}
### END TISSUEWIDE REGULATORS ###


### BEGIN TISSUESPECIFIC REGULATORS ###
pdf(file.path(outdir, "tissuespec.reg.pdf"))
act.long <- subset(act.long, tissue != "WFAT")
jmotifs.lst <- list("BFAT"="MEF2.A.B.C.D..p2", "Adr"="HNF4A_NR2F1.2.p2", "Mus"="SPIB.p2")
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

# DO OVERLAY
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

dev.off()
### END TISSUESPECIFIC REGULATORS ###

# print output
pdf(file.path(outdir, paste0(plot.i, ".regulation.pdf")))
plot.i <- plot.i + 1
multiplot(eigens.act$u.plot, eigens.act$v.plot, layout = jlayout)
dev.off()

# Regulation by DHS analysis ----------------------------------------------

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

gene.types <- c("jgenes", "jgenes.flat.filt")

jsub <- subset(df.out.lst.merged, gene.type %in% gene.types)
newlab <- hash(as.character(subset(jsub, gene.type %in% gene.types)$tissue), paste0(subset(jsub, gene.type %in% gene.types)$tissue, "\nN=", subset(jsub, gene.type %in% gene.types)$total.genes))
jsub$Tissue <- sapply(as.character(jsub$tissue), function(tiss) newlab[[tiss]])
jsub <- OrderDecreasing(jsub, "Tissue", "frac.n.spec.by.gene")
plots.dhs.sites <- ggplot(jsub, aes(x = Tissue, y = frac.n.spec.by.gene, fill = gene.type)) + geom_bar(stat = "identity", position = "dodge") + xlab("") + ylab("Fraction of genes with tissue-specific DHS within 10kb") + 
  theme_bw(18) + theme(aspect.ratio = 1, legend.position = "bottom") + scale_fill_discrete(name="Rhythmicity of Gene",
                                                                                           labels=c("Rhythmic", "Flat"))

### BEGIN MOTIF ANALYSIS ###
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
### END MOTIF ANALYSIS ###


### BEGIN WT-KO ANALYSIS ###
load("Robjs/wtko.dat.wtko.hog.commongenes.Robj", v=T)
load("Robjs/wtko.fits.all.long.commongenes.Robj", v=T)
load("Robjs/wtko.dat.complex.wtko.hog.Robj", v=T)
genes.liv <- as.character(subset(fits.best, model == "Liver")$gene)
genes.tw <- as.character(subset(fits.best, n.rhyth >= 8)$gene)
fits.sub.liv <- subset(fits.all.long.wtkohog, gene %in% genes.liv & model != "")
fits.sub.tw <- subset(fits.all.long.wtkohog, gene %in% genes.tw & model != "")

genes.liv.wtliv <- as.character(subset(fits.sub.liv, model %in% c("WT,Liver", "WT;Liver", "Liver", "WT"))$gene)
genes.liv.tw <- as.character(subset(fits.sub.tw, model %in% c("WT,Liver", "WT;Liver", "Liver", "WT"))$gene)
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
s.liv.tw <- SvdOnComplex(subset(dat.complex.wtko.hog, gene %in% genes.liv.tw), value.var = "exprs.transformed")
eigens.liv <- GetEigens(s.liv.wtliv, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2, peak.to.trough = TRUE)
eigens.tw <- GetEigens(s.liv.tw, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# multiplot(eigens.tw$v.plot, eigens.tw$u.plot, layout = jlayout)
### END WT-KO ANALYSIS ###

pdf(file.path(outdir, paste0(plot.i, ".dhs.pdf")))
plot.i <- plot.i + 1
  print(plots.dhs.sites)
  PlotLdaOut(out, jcex = 0.6, jtitle = "")
  plot(out.3$discrim[, 1], out.3$discrim[, 2], pch = ".", 
       xlab = "Liver-specific DHS vs Nonliver-specific DHS", ylab = "Liver-specific rhythmic DHS vs Others", 
       main="Motifs separating btwn tissues (x-axis) and rhythmicity (y-axis)")
  text(out.3$discrim[, 1], out.3$discrim[, 2], names(out.3$x), cex = jsize)
  abline(v = 0); abline(h = 0)
  multiplot(eigens.liv$v.plot, eigens.liv$u.plot + xlim(0, 2), layout = jlayout)
  multiplot(eigens.tw$v.plot, eigens.tw$u.plot, layout = jlayout)
dev.off()

# Liver 4CSeq -------------------------------------------------------------


