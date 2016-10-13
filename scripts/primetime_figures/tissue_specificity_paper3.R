# 2016-08-16
# Figures for paper: include hogenesch, liver kidney WTKO, nuclear proteomics, 4c-seq
# Jake Yeung


rm(list=ls())
start <- Sys.time()

library(ggplot2)
library(PMA)
# detach("package:plyr", unload=TRUE)
library(dplyr)
library(parallel)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
# source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/FourierFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/ModelStrToModel.R")
source("scripts/functions/ProteomicsFunctions.R")
source("scripts/functions/CosSineFunctions.R")

dotsize <- 6
zscore.min <- 1.25
omega <- 2 * pi / 24
n <- 4  # amplitude scaling 
# Functions ---------------------------------------------------------------



# Inits -------------------------------------------------------------------

remove.kidney.outliers <- TRUE
remove.wfat <- TRUE
plot.i <- 1
tissue.order <- c('Liver','BFAT','Kidney','Lung','Adr','Mus','Heart','Aorta','Hypo','Cere','BS')

plot.dir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots_full_paper"
dir.create(plot.dir, showWarnings = FALSE)


tfs <- GetTFs(get.mat.only = TRUE)

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch
# model selection with g1000
load("Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj", v=T)
load("Robjs/dat.complex.fixed_rik_genes.Robj", v=T)
if (remove.wfat){
  dat.complex <- subset(dat.complex, tissue != "WFAT")
  dat.long <- subset(dat.long, tissue != "WFAT")
}

fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)


# Remove kidney outliers (optional)
if (remove.kidney.outliers){
  # Kidney_SV129 genes contain some weird outliers, remove them
  outliers <- c("Cox8b", "Ucp1", "Cidea", "Flg", "Nr4a2")
  fits.long.filt <- subset(fits.long.filt, !gene %in% outliers)
  dat.wtko <- subset(dat.wtko, !gene %in% outliers)
}
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
dat.wtko.collapsed <- CollapseTissueGeno(dat.wtko)

# load proteomics

prot.long <- LoadProteomicsData()

jgenes <- c("Arntl", "Dbp", "Nr3c1", "Hsf1", "Nr1d1", "Nr1d2", "Rora", "Rorc", "Hic1", "Nrf1", "Irf2")
pdf(file.path(plot.dir, paste0("00", ".proteomics_examples.pdf")))
for (jgene in jgenes){
  PlotProteomics(subset(prot.long, gene == jgene), jtitle = jgene)
}
dev.off()


# Plot examples -----------------------------------------------------------

jgenes <- c("Slc45a3", "Insig2", "Slc44a1", "Pik3ap1", "Jun", "Mafb", "Egr1", "Nop56", "Gck", "Lipg", "Upp2", "Loxl4", "Lpin1", "Cebpb", "Pik3r1", "Hes6")

pdf(file.path(plot.dir, paste0("0", ".gene_exprs_examples.pdf")))
for (g in jgenes){
  print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == g), split.by = "tissue", jtitle = g))
}
dev.off()

# Tissue-wide modules: hogenesch -----------------------------------------------------

# get input genes
genes.tw <- as.character(subset(fits.long, n.rhyth >= 8)$gene)
genes.tw.wtko <- as.character(subset(fits.long.filt, model %in% c("Liver_SV129,Kidney_SV129"))$gene)

# get regulators: hogenesch 
outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.hogenesch"
outmain <- file.path(outbase, paste0("promoters.tissuewide.filteramp.0.15.mat"))
indir <- file.path(outmain, "expressed_genes_deseq_int.centeredTRUE")
act.long <- LoadActivitiesLong(indir, shorten.motif.name = TRUE)
act.complex <- ProjectWithZscore(act.long, omega, n = n)
sig.motifs <- unique(as.character(subset(act.complex, zscore > zscore.min)$gene))
s.act <- SvdOnComplex(subset(act.complex, tissue != "WFAT" & gene %in% sig.motifs), value.var = "exprs.transformed")

# get regulators: liver kidney WT
indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.Liver_SV129,Kidney_SV129.g=1001/atger_with_kidney.bugfixed"
act.wtko <- LoadActivitiesLongKidneyLiver(indir, shorten.motif.name = TRUE)
act.complex.wtko <- ProjectWithZscore(act.wtko, omega, n = n)
sig.motifs <- unique(as.character(subset(act.complex, zscore > zscore.min)$gene))
s.act.wtko <- SvdOnComplex(subset(act.complex.wtko, gene %in% sig.motifs), value.var = "exprs.transformed")

# plot 
comp <- 1
pdf(file.path(plot.dir, paste0(plot.i, ".tissue_wide_genes_and_regulators.pdf")))
plot.i <- plot.i + 1

s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.tw, period = 24, comp = comp, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
print(eigens.tw$u.plot)
print(eigens.tw$v.plot)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)

eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = dotsize, label.n = 25, jtitle = "", peak.to.trough = TRUE, dot.col = "black", dotsize = 2, dotshape = 18)
print(eigens.act$u.plot)
print(eigens.act$v.plot)
multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)

s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw.wtko), value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = comp, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
print(eigens$u.plot)
print(eigens$v.plot)
multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)

eigens.act <- GetEigens(s.act.wtko, period = 24, comp = comp, adj.mag = TRUE, constant.amp = dotsize, label.n = 25, jtitle = "", peak.to.trough = TRUE, dot.col = "black", dotsize = 2, dotshape = 18)
print(eigens.act$u.plot)
print(eigens.act$v.plot)
multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
dev.off()


# Hogenesch supplementary figures -----------------------------------------

# PCA
# TODO

# F24 across periods 
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

pdf(file.path(plot.dir, paste0(plot.i, ".fourier_across_periods.pdf")))
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


# Genome-wide amplitudes: Hogenesch

# from fourier directory total_variance.noise_floor.R
linesize.gwide <- 2
pdf(file.path(plot.dir, paste0(plot.i, ".genomewide_amplitude.pdf")))
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
# make Liver and Kidney the first two
jtisses <- c("Liver", "Kidney"); jlevs <- as.character(ngenes.sum$tissue)
livkid.levels <- c(jlevs[jlevs %in% jtisses], jlevs[! jlevs %in% jtisses])
dat.fit.24.ngenes.thres$tissue <- factor(as.character(dat.fit.24.ngenes.thres$tissue), levels = livkid.levels)

cbPalette <- c(gg_color_hue(2), "#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(subset(dat.fit.24.ngenes.thres, rhyth == 24), aes(x = 2 * amp.thres, y = n.genes, colour = tissue)) + 
  geom_line(size = linesize.gwide) + 
  theme_bw(24) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank()) +
  xlab("Log2 Fold Change") + ylab("# Genes") + xlim(c(0, 5)) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000)) + 
  geom_vline(xintercept = 1, linetype = "dotted") +  
  scale_colour_manual(values = cbPalette) + 
  ggtitle(paste("pval cutoff", pval.cutoff))

# Genome-wide amplitudes: Liver WT KO 
fits.bytiss <- fits.bytiss[order(fits.bytiss$amp, decreasing = TRUE), ]
fits.bytiss <- fits.bytiss  # [t]issues and [g]enotypes separate
fits.bytiss$tiss <- sapply(as.character(fits.bytiss$tissue), function(tiss) strsplit(tiss, "_")[[1]][[1]])
fits.bytiss$geno <- sapply(as.character(fits.bytiss$tissue), function(tiss) strsplit(tiss, "_")[[1]][[2]])

amp.thres <- seq(from = 0, to = max(fits.bytiss$amp), by = 0.15)
pval.cutoff <- 0.01
fits.bytiss.ngenes.thres <- subset(fits.bytiss, pval < pval.cutoff) %>%
  group_by(tiss, geno) %>%
  do(NGenesByAmp.long(., amp.thres))
fits.bytiss.ngenes.thres$rhyth <- as.factor(24)

# order by total genes
ngenes.sum.livkidWTKO <- fits.bytiss.ngenes.thres %>%
  group_by(tissue) %>%
  summarise(total = sum(n.genes)) %>%
  arrange(desc(total))
fits.bytiss.ngenes.thres$tissue <- factor(as.character(fits.bytiss.ngenes.thres$tissue), levels = ngenes.sum.livkidWTKO$tissue)
# convert "SV129" to "WT"
fits.bytiss.ngenes.thres$geno <- gsub("SV129", "WT", x = as.character(fits.bytiss.ngenes.thres$geno))
fits.bytiss.ngenes.thres$geno <- factor(as.character(fits.bytiss.ngenes.thres$geno), levels = c("WT", "BmalKO"))

ggplot(subset(fits.bytiss.ngenes.thres, rhyth == 24), aes(x = 2 * amp.thres, y = n.genes, colour = tiss, linetype = geno)) + 
  scale_linetype_manual(values = c("solid", "dashed")) + 
  geom_line(size = linesize.gwide) + 
  theme_bw(24) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank()) +
  xlab("Log2 Fold Change") + ylab("# Genes") + xlim(c(0, 5)) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000)) + 
  geom_vline(xintercept = 1, linetype = "dotted") + 
  ggtitle(paste("pval cutoff", pval.cutoff))

dev.off()


# Examples to motivate 

pdf(file.path(plot.dir, paste0(plot.i, ".gene_exprs_examps.pdf")))
plot.i <- plot.i + 1

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

library(hash)

jgenes <- c("Dbp", "Ndrg1", "Pi4k2a", "Slc44a1")
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


pdf(file.path(plot.dir, paste0(plot.i, ".model_selection_summary.pdf")))
plot.i <- plot.i + 1
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")

fits.best$n.rhyth.fac <- as.factor(sapply(as.numeric(fits.best$n.rhyth), function(n) NrhythToStr(n)))
ggplot(subset(fits.best, n.rhyth != 0), aes(x = as.factor(n.rhyth.fac), y = 2 * amp.avg)) + geom_boxplot() + xlab("# rhythmic tissues") + ylab("Log fold change")  + theme_bw(16) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# summarize number of rhyth tiss
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
  } else if (n <= 3 & n >= 2){
    return("2-3")
  } else if (n <= 7 & n >- 4){
    return("4-7")
  } else {
    print(n)
    warning("Didnt fit any if statements")
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
  xlab("Avg Log2 Fold Change") + ylab("# Genes") + xlim(c(0.15, 6)) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000)) + 
  geom_vline(xintercept = 2.8, linetype = "dotted") + 
  scale_colour_brewer(palette = "Spectral")
dev.off()


# Tissue-wide modules: liver kidney ---------------------------------------

# Find HIC1 module 

# Won't get the CIRBP module

pdf(file.path(plot.dir, paste0(plot.i, ".tissuewide_modules_liver_kidney_wtko.pdf")))
plot.i <- plot.i + 1

# BEGIN: Penalized LDA to separate between candidate TFs 
load("Robjs/liver_kidney_atger_nestle/systems_clockdriven_tissuewide_genes.Robj", v=T)


jlambda <- 0.035  # liv only
plda.out <- PenalizedLDA(M, M.labs, lambda = jlambda, K = 1, standardized = FALSE)

# plot pretty
vec.length <- sqrt(plda.out$discrim[, 1]^2)

jsize.cutoff <- 0
jsize.pairs.cut <- sapply(vec.length, function(jsize){
  if (jsize > jsize.cutoff){
    return(jsize)
  } else {
    return(0)
  }
})

labels <- names(plda.out$x)
labels.cut <- mapply(function(jlab, jsize){
  if (jsize <= 0){
    return("")
  } else {
    return(jlab)
  }
}, labels, jsize.pairs.cut)

dat.plot <- data.frame(discrim = plda.out$discrim[, 1],
                       motif = labels.cut,
                       vec.length = vec.length,
                       vec.length.cut = jsize.pairs.cut) %>% 
  mutate(discrim.floor = ifelse(discrim > 0, "Systemic", "Clock")) %>%
  arrange(discrim) %>%
  mutate(Index = seq(length(discrim)))

dat.labs <- subset(dat.plot, vec.length.cut > 0)

m <- ggplot(dat.plot, aes(x = discrim.floor, y = discrim, label = motif)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(size = 7) + 
  theme_bw(24) + 
  theme(aspect.ratio = 0.33, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("") + ylab("Motif loadings") + 
  theme(aspect.ratio = 1)
print(m)

m.index <- ggplot(dat.plot, aes(x = Index, y = discrim, label = motif)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(size = 7) + 
  theme_bw(24) + 
  theme(aspect.ratio = 0.33, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("") + ylab("Motif loadings") + 
  theme(aspect.ratio = 1, 
        axis.ticks.x = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
print(m.index)

gene.plot <- data.frame(proj = plda.out$xproj, 
                        gene = rownames(plda.out$x),
                        jlabel = ifelse(plda.out$y == 1, "Clock", "Systems"))

mm <- ggplot(gene.plot, aes(y = proj, x = as.factor(jlabel), label = gene)) + 
  geom_boxplot() +
  geom_text() + 
  theme_bw(24) + 
  xlab("") + 
  ylab("Projection") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(mm)


# END: PLDA 


# jmod1 <- "many_modules_minrhyth.4.exclude_clockdriven_model"
jmod1 <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Liver_BmalKO.Kidney_SV129,Kidney_BmalKO"
# jmod2 <- "Liver_SV129,Kidney_SV129"
jmod2 <- "Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129"
jmods <- c(jmod1, jmod2)
for (jmod in jmods){
  # outmain <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.Kidney_SV129,Kidney_BmalKO.g=1001"
  outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod, ".g=1001")
  indir <- file.path(outmain, "atger_with_kidney.bugfixed")
  act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)
  act.complex <- ProjectWithZscore(act.long, omega, n)
  sig.motifs <- unique(as.character(subset(act.complex, zscore > zscore.min)$gene))
  s.act <- SvdOnComplex(subset(act.complex, gene %in% sig.motifs), value.var = "exprs.transformed")
  
  max.labs <- 20
  jtitle <- ""
  comp <- 1
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = dotsize, 
                          label.n = max.labs, jtitle = jtitle, peak.to.trough = TRUE, label.gene = c("bHLH_family", "RORA", "SRF", "HSF1.2"), jsize = 16, 
                          dot.col = "black", dotsize = 2, dotshape = 18)
  print(eigens.act$u.plot + ylab("ZT") + ggtitle(""))
  print(eigens.act$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))
  
  # plot gene exprs module
  if (grepl("^many_modules_minrhyth", jmod)){
    n.rhyth.min <- as.numeric(strsplit(jmod, "\\.")[[1]][[2]])
    genes.tw <- as.character(subset(fits.long.filt, n.rhyth >= n.rhyth.min)$gene)
  } else {
    jmod.long <- ModelStrToModel(jmod)
    genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% jmod.long)$gene)
  }
  s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, jsize = 16)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  print(eigens$u.plot + ylab("ZT") + ggtitle(""))
  print(eigens$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))
  # plot clock rhythms
  if (jmod == jmod2){
    jmotifs.clock <- c("RORA", "bHLH_family", "HIC1", "DBP", "NFIL3")
  } else if (jmod == jmod1){
    jmotifs.clock <- c("HSF1.2", "SRF", "HIF1A", "PITX1..3", "HNF4A_NR2F1.2", "EGR1..3")
  }
  # plot example motifs
  for (m in jmotifs.clock){
    print(PlotActivitiesWithSE(subset(act.long, gene == m), jtitle = "", showSE = TRUE, jxlab = "ZT", jsize = 16) + ggtitle(m) + theme(legend.position = "none"))
  }
}

dev.off()



# Circadian flucutations in cellular heterogeneity ------------------------


exprs.thres <- 5
show.top.n <- 20
# ggplot(dat.mean.rnaseq, aes(x = exprs.mean)) + geom_density() + facet_wrap(~tissue)  # to find threshold

filt.tiss <- "WFAT"

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long.hog.genomewide <- LoadActivitiesLong(indir, shorten.motif.name = TRUE)

dat.mean.rnaseq <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))

pdf(file.path(plot.dir, paste0(plot.i, ".tissue_specific_hogenesch.pdf")))
plot.i <- plot.i + 1
jtiss <- "BFAT"
BFAT.genes <- as.character(subset(fits.long, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, BFAT.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "array"), BFAT.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, BFAT.genes, jtiss, avg.method = "mean")

jtiss <- "Mus"
Mus.genes <- as.character(subset(fits.long, model == jtiss)$gene)
m <- PlotOverlayTimeSeries(dat.long, Mus.genes, tissues = jtiss, jscale = T, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m)
m2 <- PlotOverlayTimeSeries(subset(dat.long, experiment == "array"), Mus.genes, tissues = jtiss, jalpha = 0.05, jtitle = paste0(jtiss, "-specific rhythmic genes"))
print(m2)
PlotMeanExprsOfModel(dat.mean.rnaseq, Mus.genes, jtiss, avg.method = "median")

act.long <- subset(act.long, !tissue %in% filt.tiss)
jmotifs.lst <- list("BFAT"="MEF2.A.B.C.D.", "Mus"="SPIB")
for (jtiss in names(jmotifs.lst)){
  jmotif <- jmotifs.lst[[jtiss]]
  for (jexp in c("array", "rnaseq")){
    m <- PlotActivitiesWithSE(subset(act.long.hog.genomewide, gene == jmotif & tissue == jtiss & experiment == jexp)) + theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
    m1 <- PlotActivitiesWithSE(subset(act.long.hog.genomewide, gene == jmotif & experiment == jexp)) + theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
    print(m)
    print(m1)
  }
}

# Do GO enrichment

# for (jtiss in c("BFAT", "Mus")){
# }

bfatmus.lst <- list("BFAT", "Mus")
mclapply(bfatmus.lst, function(jtiss){
  source("scripts/functions/AnalyzeGeneEnrichment.R")
  genes.bg <- as.character(subset(dat.mean.rnaseq, exprs.mean > exprs.thres & tissue == jtiss)$gene)
  genes.fg <- as.character(subset(fits.long, model == jtiss)$gene)
  enrichment <- AnalyzeGeneEnrichment(genes.bg, genes.fg)
  enrichment$minuslogpval <- -log10(as.numeric(enrichment$classicFisher))
  enrichment <- OrderDecreasing(enrichment, jfactor = "Term", jval = "minuslogpval")
  enrichment <- enrichment[1:show.top.n, ]
  m1 <- ggplot(enrichment, aes(x = Term, y = minuslogpval)) + geom_bar(stat = "identity") + 
    ylab("-log10(P-value), Fisher's exact test") + 
    xlab("") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)
  return(m1)
}, mc.cores = 2)
dev.off()



# Summarize modules -------------------------------------------------------

# filter top models
# fits.count <- subset(fits.long.filt, method == jmeth & model != "") %>% group_by(model) %>% summarise(model.count = length(model))
# fits.count <- fits.count[order(fits.count$model.count, decreasing = TRUE), ]
# fits.count <- fits.count[1:10, ]
# top.models <- as.character(fits.count$model)

# filter by hard coding
top.models <- c("Liver_SV129,Liver_BmalKO", "Liver_SV129", "Liver_BmalKO", "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO", "Liver_SV129,Kidney_SV129", "Kidney_SV129")

fits.sub <- subset(fits.long.filt, model %in% top.models)
# reorder factor by fits.count
fits.sub$model <- factor(as.character(fits.sub$model), levels = top.models)


amp.thres <- seq(from = 0, to = max(fits.bytiss$amp), by = 0.15)
fits.counts.by.amp <- fits.sub %>%
  group_by(model) %>%
  do(NGenesByAmp.long(., amp.thres, labelid = "model", varid = "amp.avg", outlabel = "model"))

pdf(file.path(plot.dir, paste0(plot.i, ".summarize_livkid_modules.pdf")))
plot.i <- plot.i + 1
mm <- ggplot(fits.counts.by.amp, aes(x = 2 * amp.thres, y = n.genes, group = model, colour = as.factor(model))) + geom_line() +
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
print(mm)
print(mm + theme(legend.position = "none"))
dev.off()



# Single TFs underlie clock-independent liver-specific diurnal gen --------

dat.mean.wtko <- dat.wtko.collapsed %>%
  group_by(gene, tissue, geno) %>%
  summarise(exprs.mean = mean(exprs))

# jmod <- "Liver_SV129,Liver_BmalKO"

jmods <- c("Liver_SV129,Liver_BmalKO", "Liver_BmalKO", "Kidney_SV129", "Liver_SV129")
for (jmod in jmods){
  jtiss <- strsplit(jmod, "_")[[1]][[1]]
  print(paste("Finding regulators for:", jmod))
  genes.mod <- unique(as.character(subset(fits.long.filt, model == jmod)$gene))
  
  maraoutdir <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod, ".g=1001/atger_with_kidney.bugfixed")
  act.s <- LoadActivitiesLongKidneyLiver(maraoutdir, shorten.motif.name = TRUE)
  act.s.complex <- ProjectWithZscore(act.s, omega, n)
  sig.motifs <- unique(as.character(subset(act.s.complex, zscore > zscore.min)$gene))
  s.act <- SvdOnComplex(subset(act.s.complex, gene %in% sig.motifs), value.var = "exprs.transformed")
  
  
  jmod.label <- gsub(pattern = ",", replacement = "-", jmod)
  pdf(file.path(plot.dir, paste0(plot.i, ".clock_independent_", jmod.label, "_regulators.pdf")))
  

  # increment plot.i after the for loop
  
  s <- SvdOnComplex(subset(dat.freq, gene %in% genes.mod), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = comp, label.n = 20, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, label.gene = c("Mafb", "Egr1", "Creb3", "Elf2"))
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = dotsize, label.n = 20, jtitle = "", peak.to.trough = TRUE, dot.col = "black", dotsize = 2, dotshape = 18, label.gene = c("ELF1.2.4"))
  
  # plot genes and regulators
  print(eigens$u.plot + ylab("ZT") + ggtitle(""))
  print(eigens$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
  print(eigens.act$u.plot + ylab("ZT") + ggtitle(""))
  print(eigens.act$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))
  multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
  # show mean exprs across tissues
  print(PlotMeanExprsOfModel(dat.mean.wtko, genes = genes.mod, jmodel = jmod, sorted = TRUE, avg.method = "mean"))
  
  # Print genes that match model
  jmotifs <- names(head(eigens.act$eigensamp[order(abs(eigens.act$eigensamp), decreasing = TRUE)], n = 20))
  for (jmotif in jmotifs){
    genes.all <- unlist(sapply(jmotif, GetGenesFromMotifs, tfs))
    genes.that.fit <- as.character(subset(fits.long.filt, gene %in% genes.all & model %in% jmod & method == jmeth)$gene)
    if (length(genes.that.fit) > 0){
      for (gene.hit in genes.that.fit){
        if (jmod %in% c("Liver_SV129,Liver_BmalKO", "Liver_SV129", "Liver_BmalKO")){
          gene.prot <- gene.hit
          jprot.long <- prot.long
        } else {
          gene.prot <- ""
          jprot.long <- NA
        }
        print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == gene.hit), jtitle = gene.hit))
        print(PlotActivitiesWithSE(subset(act.s, gene == jmotif), jtitle = jmotif) + theme_bw())
        print(PlotmRNAActivityProtein(dat.wtko, act.s, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22))
        print(PlotmRNAActivityProtein(dat.wtko, act.s, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = "both", dotsize = 2, themesize = 14) + theme(strip.text = element_blank()))
      }
    }
  }

  dev.off()  
}
plot.i <- plot.i + 1




# Go analysis on liver kidney WT KO modules -------------------------------
# Optional: go_terms_analysis_with_phase_amplitude.R provides a cooler picture

# # Do GO enrichment
# pdf(file.path(plot.dir, paste0(plot.i, ".go_analysis.pdf")))
# 
# jtiss.lst <- list(ModelStrToModel(jmod1),
#                   ModelStrToModel(jmod2),
#                     "Liver_SV129", 
#                     "Liver_SV129,Liver_BmalKO", 
#                     "Kidney_SV129", 
#                     "Kidney_SV129,Kidney_BmalKO", 
#                     "Liver_BmalKO")
# 
# # for (jtiss in c("Liver_SV129,", "Liver_SV129", "Liver_SV129,Liver_BmalKO", "Kidney_SV129", "Kidney_SV129,Kidney_BmalKO", "Liver_BmalKO")){
# # }
# 
# 
# 
# 
# jtiss.lst <- jtiss.lst[1]
# mclapply(jtiss.lst, function(jtiss){
#   source("scripts/functions/AnalyzeGeneEnrichment.R")
#   source("/home/yeung/projects/tissue-specificity/scripts/functions/ListFunctions.R")
#   plot.lst <- expandingList()
# # lapply(jtiss.lst, function(jtiss){
#   print(jtiss)
#   genes.bg <- as.character(subset(fits.long.filt)$gene)
#   genes.fg <- as.character(subset(fits.long.filt, model %in% jtiss)$gene)
#   # genes.fg <- genes.fg[! genes.fg %in% go.genes]  # remove genes from one go.term to see if we can get others
#   # genes.filt <- c("Trdmt1", "Trmt5", "Mettl1", "Nsun2", "Trtm61a")  # tRNA methylation??
#   enrichment <- AnalyzeGeneEnrichment(genes.bg, genes.fg, FDR.cutoff = 0.5, return.GOdata = TRUE)
#   enrichment$minuslogpval <- -log10(as.numeric(enrichment$classicFisher))
#   enrichment <- OrderDecreasing(enrichment, jfactor = "Term", jval = "minuslogpval")
#   show.top.n.min <- min(nrow(enrichment), show.top.n)
#   if (show.top.n.min == 0) return(NULL)
#   enrichment <- enrichment[1:show.top.n.min, ]   # prevent taking more than you have enrichment
#   m1 <- ggplot(enrichment, aes(x = Term, y = minuslogpval)) + geom_bar(stat = "identity") +
#     ylab("-log10(P-value), Fisher's exact test") +
#     xlab("") +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#     ggtitle(jtiss)
#   plot.lst$add(m1)
#   # plot genes in enrichment
#   i <- 1
#   max.genes <- 40
#   # print(show.top.n.min)
#   # print(enrichment)
#   for (i in seq(show.top.n.min)){
#     go.genes <- enrichment$genes[[i]]
#     go.term <- enrichment$Term[[i]]
#     if (is.na(go.genes)) next
#     show.n.genes <- min(length(go.genes), max.genes)
#     s <- SvdOnComplex(subset(dat.freq, gene %in% go.genes), value.var = "exprs.transformed")
#     eigens <- GetEigens(s, period = 24, comp = comp, label.n = show.n.genes, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, label.gene = c("Mafb", "Egr1", "Creb3"))
#     plot.lst$add(eigens$u.plot + ggtitle(go.term))
#     # print(eigens$v.plot + ggtitle(go.term))
#   }
#   return(plot.lst$as.list())
# # })
# }, mc.cores = length(jtiss.lst))
# dev.off()
# plot.i <- plot.i + 1


# Cooperative TFs underlie clock-dependent tissue-specific diurnal --------

# show enrichment of DHS peaks
wtmodulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.1000.tissue.Liver.module.Liver_SV129random.flat.TRUE.Robj"
load(wtmodulef, v=T)
df.out.lst.merged.liverWT <- df.out.lst.merged
df.out.lst.merged.liverWT$gene.type[1:3] <- c("Liver_SV129", "Flat_SV129", "Flat.filt_SV129")

wtko.modulef <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks/n.tiss.spec.df.out.lst.rand.1000.tissue.Liver.module.Liver_SV129,Liver_BmalKOrandom.flat.TRUE.Robj"
load(wtko.modulef, v=T)
df.out.lst.merged.liverWTKO <- df.out.lst.merged
df.out.lst.merged.liverWTKO$gene.type[1:3] <- c("Liver_SV129.Liver_BmalKO", "Flat_SV129BmalKO", "Flat.filt_SV129BmalKO")
df.out.lst.merged <- rbind(df.out.lst.merged.liverWT, df.out.lst.merged.liverWTKO)
df.out.lst.merged$xlabs <- make.names(df.out.lst.merged$gene.type, unique = TRUE)
df.out.lst.merged <- df.out.lst.merged[order(df.out.lst.merged$gene.type, decreasing = FALSE), ]

df.out.lst.bg <- df.out.lst.merged[grepl("^Random", df.out.lst.merged$xlabs), ]
df.out.lst.bg.meanvar <- df.out.lst.bg %>%
  group_by(gene.type) %>%
  summarise(mean.frac = mean(frac.n.spec.by.gene),
            sd.frac = sd(frac.n.spec.by.gene))
df.out.lst.fg.meanvar <- df.out.lst.merged[!grepl("^Random", df.out.lst.merged$xlabs), ] %>%
  group_by(gene.type) %>%
  summarise(mean.frac = mean(frac.n.spec.by.gene),
            sd.frac = sd(frac.n.spec.by.gene))
df.out.lst.meanvar <- rbind(df.out.lst.fg.meanvar[!grepl("^Flat", df.out.lst.fg.meanvar$gene.type), ], df.out.lst.bg.meanvar)

barwidth <- 0.8
limits <- aes(ymax = mean.frac + sd.frac, ymin=mean.frac - sd.frac)
m.bar <- ggplot(df.out.lst.meanvar, aes(x = gene.type, y = mean.frac)) + geom_bar(stat = "identity", width = barwidth) + theme_bw() + 
  geom_errorbar(limits, width = barwidth / 2) + theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  xlab("") + ylab("Fraction of Genes with Liver-specific DHS peaks")

jmod <- "Liver_SV129,Liver_BmalKO"
load("Robjs/penalized_lda_robjs/mat.pmd.RORA_bHLH_SRF.40000.g1001.Robj", v=T)
pldarobj <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/plda_robjs/2D.posterior.multigene.distfilt.morenonliv.bugfixed.liverWTKO.40000.cutoff.3.cutofflow0.method.g=1001.Robj"
load(pldarobj)

mat.fgbg.lab.lst.3 <- SetUpMatForLda(mat.fg, mat.bgnonliver, mat.bg, has.peaks = TRUE)
mat.fgbg.3 <- mat.fgbg.lab.lst.3$mat.fgbg; labels3 <- mat.fgbg.lab.lst.3$labels
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), RemoveP2Name)
colnames(mat.fgbg.3) <- sapply(colnames(mat.fgbg.3), function(cname){
  return(RemoveCommasBraces(cname))
}, USE.NAMES = FALSE)


rhyth.motifs <- sapply(GetTopMotifs("rhythmic"), RemoveP2Name, USE.NAMES = FALSE)
rhyth.motifs <- c(rhyth.motifs, c("SRF"))
# add DBP
# rhyth.motifs <- c(rhyth.motifs, c("DBP"))
rhyth.motifs <- rhyth.motifs[which( ! rhyth.motifs %in% c("ATF6"))]
tissue.motifs <- sapply(GetTopMotifs("tissue"), RemoveP2Name, USE.NAMES = FALSE)
tissue.motifs <- c(tissue.motifs, c("ATF5_CREB3", "ATF6"))
tissue.motifs <- tissue.motifs[which( ! tissue.motifs %in% c("SRF"))]

# remove tissue motifs in rhyth
rhyth.motifs <- rhyth.motifs[which(!rhyth.motifs %in% intersect(rhyth.motifs, tissue.motifs))]

# cross prods
mat.rhyth3 <- subset(mat.fgbg.3, select = intersect(rhyth.motifs, colnames(mat.fgbg.3)))
mat.tiss3 <- subset(mat.fgbg.3, select = intersect(tissue.motifs, colnames(mat.fgbg.3)))
mat.rhythtiss3 <- CrossProductTwoSets(mat.rhyth3, mat.tiss3)

mat.fgbg.cross.rhythtiss3 <- cbind(mat.fgbg.3, mat.rhythtiss3)
# remove columns with 0 variance
mat.fgbg.cross.rhythtiss3[which(colSums(mat.fgbg.cross.rhythtiss3) == 0)] <- list(NULL)

jlambda <- 0.035  # liv only
out.cross.rhythtiss3 <- PenalizedLDA(mat.fgbg.cross.rhythtiss3, labels3, lambda = jlambda, K = 2, standardized = FALSE)

# plot pretty
vec.length <- sqrt(out.cross.rhythtiss3$discrim[, 1]^2 + out.cross.rhythtiss3$discrim[, 2]^2)

jsize.cutoff <- 0.1
jsize.pairs.cut <- sapply(vec.length, function(jsize){
  if (jsize > jsize.cutoff){
    return(jsize)
  } else {
    return(0)
  }
})

labels <- names(out.cross.rhythtiss3$x)
labels.cut <- mapply(function(jlab, jsize){
  if (jsize <= 0){
    return("")
  } else {
    return(jlab)
  }
}, labels, jsize.pairs.cut)

dat.plot <- data.frame(x = out.cross.rhythtiss3$discrim[, 1],
                       y = out.cross.rhythtiss3$discrim[, 2],
                       motif = labels.cut,
                       vec.length = vec.length,
                       vec.length.cut = jsize.pairs.cut)
dat.labs <- subset(dat.plot, vec.length.cut > 0)

m <- ggplot(dat.plot, aes(x = x, y = y)) + 
  geom_point(size = 0.01) + 
  geom_text_repel(data = dat.labs, aes(x = x, y = y, label = motif), size = 1.5) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_bw() + 
  theme(aspect.ratio = 0.33, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Motif loadings separating liver and kidney DHS peaks") + ylab("Motif loadings separating rhythmic and flat DHS peaks")


# MARA on subset of motifs

descrip <- "sep_liv_rhyth.RORA_bHLH_SRF.fixsignbugfixed"
K <- 3

outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/mara_out.", descrip, ".K_", K)
inmain <- outmain
indir <- file.path(inmain, "atger_with_kidney.bugfixed")
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/many_motifs_decorrelated.nofiltergenes.noabs/atger_with_kidney.bugfixed"
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)
act.long$pc <- sapply(as.character(act.long$gene), function(g) as.numeric(strsplit(g, "_")[[1]][[3]]), USE.NAMES = FALSE)
motif.decor <- sapply(act.long$pc, function(p) PcToMotif(mat.pmd$v, p), USE.NAMES = FALSE)
act.long$gene <- motif.decor

act.complex <- act.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))


s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# jtitle <- gsub(pattern = "\\.", replacement = "\n", basename(indirmain))

max.labs <- 2
jtitle <- ""
eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = dotsize, 
                        label.n = K, jtitle = jtitle, peak.to.trough = TRUE, dot.col = "black", 
                        dotsize = 2, dotshape = 18, xlab = "Activity (arbitrary units)")


## Do motif cooperativity using contingency
K <- 300
prefix <- "Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K."
inf <- paste0(prefix, K, ".Robj")
load(inf, v=T)

jmotifs <- c("FOXA2", "ONECUT1.2", "CUX2")
grepstr <- paste(jmotifs, collapse = "|")
(fits.sub <- subset(fits, grepl(grepstr, pair)) %>% arrange(pval))

fits.sub$has.ROR <- sapply(fits.sub$pair, function(p) grepl("RORA", p))

fits.sub$pair <- factor(fits.sub$pair, levels = fits.sub$pair)
m.loglin <- ggplot(subset(fits.sub, pval < 0.05), aes(x = pair, y = -log10(pval), fill = has.ROR)) + geom_bar(stat = "identity") +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Top motif pairs between\nFOXA2, ONECUT, or CUX2") + 
  xlab("")


pdf(file.path(plot.dir, paste0(plot.i, ".clock_dependent_liver_genes_cooperative_action.pdf")))
plot.i <- plot.i + 1
print(m.bar)
print(m)
print(eigens.act$u.plot)
print(m.loglin)
multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
dev.off()


# 4CSeq reveals tissue-specific promoter-enhancer contacts ----------------

projdir <- "/home/yeung/projects/4c_seq"
load(file.path(projdir, "Robjs/counts.long.delt.lk.LR.Robj"), v=T)
load(file.path(projdir, "Robjs/counts.hits.ignore5frags.normInf.with_assay_signal.Robj"), v=T)

source(file.path(projdir, "scripts/functions/GetBaitLocations.R"))
source(file.path(projdir, "scripts/functions/PlotFunctions.R"))
source(file.path(projdir, "scripts/functions/TrackHubFunctions.R"))
source(file.path(projdir, "scripts/functions/AssayFunctions.R"))
source(file.path(projdir, "scripts/functions/FindPeaks.R"))
source(file.path(projdir, "scripts/functions/MergeCounts.R"))

bait.locs <- GetBaitLocations()

zt <- "ZT20"
zt08 <- "ZT08"
jdist <- 250000

jbaits <- c("Mreg", "GYS2", "Pi3kap1", "Slc44a1", "Slc45a3short", "Slc45a3long")
baits.newnames <- c("Mreg", "Gys2", "Pik3ap1", "Slc44a1", "Slc45a3 (short)", "Slc45a3 (long)")
baits.genename <- c("Mreg", "Gys2", "Pik3ap1", "Slc44a1", "Slc45a3", "Slc45a3")
# jbaits <- c("Mreg", "Pi3kap1", "Slc44a1", "Slc45a3short", "Slc45a3long")
# baits.newnames <- c("Mreg", "Pik3ap1", "Slc44a1", "Slc45a3", "Slc45a3")
baits.hash <- hash(jbaits, baits.newnames)
baits.genename.hash <- hash(jbaits, baits.genename)

remove.bait <- c("Gys2", "Slc45a3 (long)")
# rename for volcano plots
counts.delt.lk.sub$bait <- factor(sapply(as.character(counts.delt.lk.sub$bait), function(b) baits.hash[[b]]), levels = unique(baits.newnames))

pdf(file.path(plot.dir, paste0(plot.i, ".4cseq_tissue_plots.nogys2.pdf")))
plot.i <- plot.i + 1

for (jbait in jbaits){
  jbait.new <- baits.hash[[jbait]]
  jbait.genename.new <- baits.genename.hash[[jbait]]
  print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == jbait.genename.new), jtitle = jbait.genename.new))
  
  jsub.sig <- subset(counts.long, genotype == "WT" & time == zt)
  jsub.sig08 <- subset(counts.long, genotype == "WT" & time == zt08)
  Signal <- PlotSignalLivVsKidLR(jbait, jsub.sig, pseudo.low = 500, jtitle = paste(jbait.new, zt), mindist = jdist, do.facet = FALSE)
  Signal08 <- PlotSignalLivVsKidLR(jbait, jsub.sig08, pseudo.low = 500, jtitle = paste(jbait.new, zt08), mindist = jdist, do.facet = FALSE, show.legend = FALSE)
  Zscores <- MergeCountsLivKid(jbait, counts.delt.lk, bait.locs, max.dist = jdist, show.plot = "Zscore", jtitle = "", jxlab = "", jshow.legend = FALSE, flip.y.axis = TRUE)
  Pvalues <- MergeCountsLivKid(jbait, counts.delt.lk, bait.locs, max.dist = jdist, show.plot = "Pvalue", jtitle = "", jshow.legend = TRUE, flip.y.axis = TRUE)
  PlotQuad(Signal, Signal08, Zscores, Pvalues, n.boxes = 2)
  PlotTriple(Signal, Signal08, Zscores + theme(aspect.ratio = 0.25), n.boxes = 2)
}

ggplot(subset(counts.delt.lk.sub, bait == "Mreg"), aes(x = A.delta / log10(2), y = -log10(pval.row))) +
  theme_bw() + 
  geom_point(alpha = 0.5) + facet_grid(time ~ bait) + xlim(-2.4, 2.4) + 
  xlab("Log fold change (Kidney - Liver)") + 
  ylab("-log10(p-value") + 
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(counts.delt.lk.sub, bait %in% baits.newnames & ! bait %in% remove.bait), aes(x = A.delta / log10(2), y = -log10(pval.row))) +
  theme_bw() + 
  geom_point(alpha = 0.5) + facet_grid(time ~ bait) + xlim(-2.4, 2.4) + 
  xlab("Log fold change (Kidney - Liver)") + 
  ylab("-log10(p-value") + 
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(counts.delt.lk.sub, bait %in% baits.newnames), aes(x = A.delta / log10(2), y = -log10(pval.row))) +
  theme_bw() + 
  geom_point(alpha = 0.5) + facet_grid(time ~ bait) + xlim(-2.4, 2.4) + 
  xlab("Log fold change (Kidney - Liver)") + 
  ylab("-log10(p-value") + 
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dev.off()



# Alternative promoter ----------------------------------------------------

# TODO

print(Sys.time() - start)