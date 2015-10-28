# 2015-10-20 
# Figures that will be part of paper to analyze this hogenesch dataset

remove.wfat <- TRUE

tissue.order <- c('Liver','BFAT','Kidney','Lung','Adr','Mus','Heart','Aorta','Hypo','Cere','BS')

setwd("/home/yeung/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots_paper"

library(ggplot2)
library(hash)
library(dplyr)
library(reshape2)
library(gplots)
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
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/OuterComplex.R")
source("scripts/functions/GetTFs.R")




# Figure 1A Present high level variation BETWEEN TISSUES and WITHIN TISSUES (time) ------------------
# scripts from pca_adjusted_microarray.label_variance.R

# Load data, log transform

log2.transform <- TRUE
load("Robjs/dat.array.adj.primetime.Robj", verbose=T)  # no -Infs

# dat <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt", 
#                            remove.negs = TRUE, fix.rik.xgene = TRUE)
# dat <- LoadArray(form = "wide")  # not adjusted, set log2.transform = FALSE
# dat <- LoadRnaSeq(handle.duplicates = FALSE)  # do RNA-Seq instead, log2.transform = TRUE

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

# Consolidate PCA into long
jtissues <- GetTissues(rownames(dat_pca$x), get_unique = FALSE)
jtimes <- as.numeric(GetTimes(rownames(dat_pca$x), get_unique = FALSE))
pca.long <- data.frame(tissue = rep(jtissues, times = ncol(dat_pca$x)),
                       time = rep(jtimes, times = ncol(dat_pca$x)),
                       pc = rep(colnames(dat_pca$x), each = nrow(dat_pca$x)),
                       loading = unlist(as.vector(dat_pca$x), use.names = FALSE))
head(pca.long)
n.samps <- length(jtissues)

# Analyze whether each PC is rhythmic or not
pca.p <- pca.long %>%
  group_by(pc, tissue) %>%
  do(GetPeriodogramFreq(.))

# Summarize by number of tissues with 12, 24, Inf, or other
pca.p.sum <- pca.p %>%
  group_by(pc) %>%
  do(SummarisePeriodogram2(., weighted = TRUE))
# subset(pca.p.sum, pc == "PC13")
# sum(subset(pca.p.sum, pc == "PC13")[1, 2:5])

pca.p.sum$pc.num <- sapply(pca.p.sum$pc, function(p) as.numeric(strsplit(as.character(p), split = "PC")[[1]][[2]]))
pca.p.sum <- pca.p.sum[order(pca.p.sum$pc.num), ]

# Calculate eigenvalues and attach to pca.psum
eigenvals <- dat_pca$sdev ^ 2 / sum(dat_pca$sdev ^ 2)
pcs <- c(paste("PC", seq(n.samps), sep = ""))
eigenvals.dic <- hash(pcs, eigenvals)

pca.p.sum$eigenvals <- sapply(as.character(pca.p.sum$pc), function(pc) eigenvals.dic[[pc]])
pca.p.sum$pc.num <- sapply(as.character(pca.p.sum$pc), function(x) as.numeric(substr(x, 3, nchar(x))))
# adjust factors for plotting
pca.p.sum$pc <- factor(as.character(pca.p.sum$pc), levels = pcs)

# Make long and plot
pca.p.sum.long <- melt(pca.p.sum, id.var = c("pc", "eigenvals", "pc.num"), value.name = "fracFourier", variable.name = "Component")

pca.p.sum.long <- pca.p.sum.long %>%
  mutate(eigenvals.frac = eigenvals * fracFourier)

# plot PCA from 1 to 12, 13 to 100
starts <- c(1, 13, 100)
ends <- c(12, 50, 288)
pca.p.sum.long$Component <- factor(pca.p.sum.long$Component, levels = c("T.inf", "T.other",  "T.12", "T.24"))
pca.p.sum.long <- pca.p.sum.long %>%
  group_by(pc) %>%
  arrange(Component)


pdf(file.path(outdir, "PCA_scree.pdf"))
for (i in seq(length(starts))){
  print(PlotComponents2(pca.p.sum.long, starts[[i]], ends[[i]]))  # PCA 1,2,13,14,15,16, .. .20 are interesting
}
dev.off()

# plot first 50 modules across tissues and time

# order pca.p.sum.long tissues to same as Fourier components: hardcode
pca.long$tissue <- factor(pca.long$tissue, levels = tissue.order)

pdf(file.path(outdir, "PCA_modules.pdf"))
pcs <- seq(50); pcs <- paste("PC", pcs, sep = "")
for (jpc in pcs){
  print(ggplot(subset(pca.long, pc == jpc), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jpc) + 
          xlab("CT") +
          ylab("PC loading") + 
          theme_bw(24) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
          scale_x_continuous(breaks=c(18, 42, 64)))
}
dev.off()

# Plot interesting screeplots and PCs into a single panel
scree1 <- PlotComponents2(pca.p.sum.long, starts[[1]], ends[[1]])
scree2 <- PlotComponents2(pca.p.sum.long, starts[[2]], ends[[2]])
jpc <- "PC1"
pc1 <- ggplot(subset(pca.long, pc == jpc), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jpc) + 
                 xlab("CT") +
                 ylab("PC loading") + 
                 theme_bw(24) + 
                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
                 scale_x_continuous(breaks=c(18, 42, 64))
jpc <- "PC13"
pc2 <- ggplot(subset(pca.long, pc == jpc), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jpc) + 
  xlab("CT") +
  ylab("PC loading") + 
  theme_bw(24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
  scale_x_continuous(breaks=c(18, 42, 64))


multiplot(scree1, scree2, pc1, pc2, cols = 2)  


# Figure 1B Temporal variation using Fourier analysis ------------------
# from fourier/total_variance.R
# load("Robjs/dat.complex.all_periods.array_only.Robj", verbose=T)
load("Robjs/dat.complex.all_T.rbinded.Robj", verbose=T)

if (remove.wfat){
  dat.complex.all_T <- subset(dat.complex.all_T, tissue != "WFAT")
}

dat.var.s <- dat.complex.all_T %>%
  group_by(period, tissue) %>%
  summarise(sum_sqr_mod = sum(Mod(exprs.transformed) ^ 2)) %>%
  mutate(period.factor = signif(period, digits = 3))

dat.var.s$period.factor <- factor(dat.var.s$period.factor, 
                                  levels = sort(unique(dat.var.s$period.factor), decreasing = TRUE))

dat.var.all <- dat.var.s %>%
  group_by(tissue) %>%
  summarise(sum_sqr_mod.total = sum(sum_sqr_mod))
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

pdf(file.path(outdir, "fourier_variance.pdf"))
# for ordering the facet_wrap plot across tissues by total variance
dat.var.s$tissue <- factor(dat.var.s$tissue,
                           levels = dat.var.all$tissue[order(dat.var.all$sum_sqr_mod.total, decreasing = TRUE)])

cbPalette <- c("#009E73", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(subset(dat.var.s), aes(x = tissue, y = sum_sqr_mod, fill = period.factor.cond)) + 
  geom_bar(stat="identity") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") + xlab("") + ylab("Spectral power") +
  scale_fill_manual(name = "Fourier component", drop=TRUE, values=cbPalette)

# for ordering the facet_wrap plot across tissues by 24h variance
dat.var.s$tissue <- factor(dat.var.s$tissue,
                           levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
ggplot(subset(dat.var.s, period.factor.cond != "Other"), aes(x = tissue, y = sum_sqr_mod, fill = period.factor.cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") + xlab("") + ylab("Spectral power") +
  scale_fill_manual(name = "Fourier component", drop=TRUE, values=cbPalette)
dev.off()

pdf(file.path(outdir, "circadian_and_ultraidian_variance.pdf"))
# plot normalized spectral power
dat.var.s1_adj$tissue <- factor(dat.var.s1_adj$tissue,
                                levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])
ggplot(subset(dat.var.s1_adj, period.factor.cond != "Other"), aes(x = tissue, y = s1_normalized, fill = period.factor.cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") + xlab("") + ylab("Normalized spectral power") +
  scale_fill_manual(name = "Fourier component", drop=TRUE, values=cbPalette)

ggplot(dat.var.s1_adj.24, aes(x = tissue, y = s1_normalized)) + geom_bar(stat = "identity") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") + xlab("") + ylab("Normalized 24h spectral power") + 
  facet_wrap(~period.factor.cond)

# plot normalized spectral power
ggplot(dat.var.s1_adj.12, aes(x = tissue, y = s1_normalized)) + geom_bar(stat = "identity") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") + xlab("") + ylab("Normalized 12h spectral power") + 
  facet_wrap(~period.factor.cond)
dev.off()

dat.var.s1_adj$tissue <- factor(dat.var.s1_adj$tissue,
                                levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])
pdf(file.path(outdir, "normalized_fourier_variance.pdf"))

m1 <- ggplot(dat.var.s1_adj, aes(x = period.factor, y = s1_normalized)) +  geom_bar(stat = "identity") + facet_wrap(~tissue) + 
  theme_bw(24) + theme(axis.text.x=element_text(size=11, angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") +
  facet_wrap(~tissue) +
  xlab("Fourier component (h)") + ylab("Normalized spectral power")
print(m1)
dev.off()

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

pdf(file.path(outdir, "fourier_across_periods.pdf"))

xscale_periods <- seq(6, 30, 2)
ggplot(subset(dat.fit.periods.sub), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.sub$period))/55) + geom_vline(xintercept=24, linetype="dotted") + 
  scale_x_continuous(breaks=xscale_periods) + 
  xlab("Period with minimum RSS (h)") + ylab("Number of genes") +
  theme_bw(24) + theme(panel.grid.major = element_blank(), 
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
  theme_bw(24) + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       axis.text.x=element_text(size=11))
print(m2)

ggplot(dat.fit.periods.sub, aes(x = period, y = ssq.residuals, colour = gene)) + geom_point() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")

# Observe Cry1 Aorta why it is period of 22
tiss <- "Aorta"; gen <- "Cry1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
w <- 2 * pi / 24
w.min <- 2 * pi / period.min
fit24 <- lm(exprs ~ 0 + experiment + sin(w * time) + cos(w * time), dat.sub)
fit.min <- lm(exprs ~ 0 + experiment + sin(w.min * time) + cos(w.min * time), dat.sub)
dat.sub.array <- subset(dat.sub, experiment == exper)
plot(dat.sub.array$time, predict(fit24, dat.sub.array), "o", col = "blue", ylim = range(dat.sub.array$exprs), 
     main = paste0(tiss, " ", gen, " ", exper, " T=24h (blue) vs T=", period.min, "h (red)"))
lines(dat.sub.array$time, predict(fit.min, dat.sub.array), "o", col = "red")
points(dat.sub.array$time, dat.sub.array$exprs, pch="*", col = "black")

# observe why Cry1 Mus has period 26.3
tiss <- "Mus"; gen <- "Cry1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Mus"; gen <- "Nr1d1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Hypo"; gen <- "Nfil3"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Hypo"; gen <- "Wee1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "BFAT"; gen <- "Myh7"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Adr"; gen <- "Arntl"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.genome_wide.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

dev.off()

pdf(file.path(outdir, "fig1_bottom.pdf"))
# Plot normalized spectral power and period with minimum RSS side by side
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(m1, m2, layout = jlayout)  
dev.off()

# Figure 2 Tissue-modules ----------------------------------------------------------

filt.tiss <- c("WFAT")
load("Robjs/fits.best.collapsed_models.Robj", verbose=T)
load("Robjs/dat.complex.fixed_rik_genes.Robj")

if (remove.wfat){
  dat.complex <- subset(dat.complex, ! tissue %in% filt.tiss)
}

pdf(file.path(outdir, "tissue_modules.tissue_wide.pdf"))

# Plot global statistics
ggplot(fits.best, aes(x = weight)) + geom_density() + facet_wrap(~n.rhyth)
ggplot(fits.best, aes(x = weight, y = amp.avg)) + geom_point(alpha = 0.5) + facet_wrap(~n.rhyth)

# Plot counts


# Plot tissue wide genes
fits.tw <- subset(fits.best, n.rhyth >= 8)
genes.tw <- as.character(fits.tw$gene)
#outobj <- PlotHeatmapNconds(fits.tw, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)

s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
eigens.tw <- GetEigens(s.tw, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)  

dev.off()

pdf(file.path(outdir, "tissue_modules.tissue_specific.pdf"))

# Plot tissue specific genes
tissues <- unique(dat.long$tissue)
for (tiss in tissues){
  fits.ts <- subset(fits.best, n.rhyth == 1 & model == tiss)
  genes.ts <- as.character(fits.ts$gene)
  print(paste("N genes for", tiss, ":", length(genes.ts)))
  if (length(genes.ts) <= length(tissues)){
    # problems with too few genes skip
    next
  }  
  s.ts <- SvdOnComplex(subset(dat.complex, gene %in% genes.ts), value.var = "exprs.transformed")
  eigens.ts <- GetEigens(s.ts, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.ts$u.plot, eigens.ts$v.plot, layout = jlayout)    
}

dev.off()

pdf(file.path(outdir, "tissue_modules.pairs_triplets.pdf"))
# Plot striking modules
# Aorta and BFAT
# fits.bfataorta <- subset(fits.best, n.rhyth >= 2 & n.rhyth <= 3)
fits.bfataorta <- subset(fits.best, n.rhyth == 2)
fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT|BFAT.*Aorta", fits.bfataorta$model), ]
# fits.bfataorta <- fits.bfataorta[grep("Aorta;BFAT|BFAT;Aorta", fits.bfataorta$model), ]
# fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT|BFAT.*Aorta", fits.bfataorta$model), ]

genes.bfataorta <- as.character(fits.bfataorta$gene)
print(paste("Genes in bfat-aorta:", length(genes.bfataorta)))

#outobj <- PlotHeatmapNconds(fits.bfataorta, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)
s.bfataorta <- SvdOnComplex(subset(dat.complex, gene %in% genes.bfataorta & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
eigens.bfataorta <- GetEigens(s.bfataorta, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.bfataorta$u.plot, eigens.bfataorta$v.plot, layout = jlayout)  

# BFAT, Aorta, Muscle
fits.bfataortamus <- subset(fits.best, n.rhyth >= 2 & n.rhyth <= 3)
fits.bfataortamus <- fits.bfataortamus[grep("Aorta.*BFAT.*Mus", fits.bfataortamus$model), ]

genes.bfataortamus <- as.character(fits.bfataortamus$gene)
print(paste("Genes in bfat-aorta-mus:", length(genes.bfataortamus)))

#outobj <- PlotHeatmapNconds(fits.bfataortamus, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)
s.bfataortamus <- SvdOnComplex(subset(dat.complex, gene %in% genes.bfataortamus & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
eigens.bfataortamus <- GetEigens(s.bfataortamus, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.bfataortamus$u.plot, eigens.bfataortamus$v.plot, layout = jlayout)  

# Liver Kidney
# fits.livkid <- subset(fits.best, n.rhyth >= 2 & n.rhyth <= 3)
fits.livkid <- subset(fits.best, n.rhyth == 2 | n.rhyth == 3)
fits.livkid <- fits.livkid[grep("Liver.*Kidney|Kidney.*Liver", fits.livkid$model), ]
genes.livkid <- as.character(fits.livkid$gene)
print(paste("Genes in liver-kidney:", length(genes.livkid)))

#outobj <- PlotHeatmapNconds(fits.livkid, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)
s.livkid <- SvdOnComplex(subset(dat.complex, gene %in% genes.livkid & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
eigens.livkid <- GetEigens(s.livkid, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.livkid$u.plot, eigens.livkid$v.plot, layout = jlayout)  

eigens.livkid <- GetEigens(s.livkid, period = 24, comp = 2, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.livkid$u.plot, eigens.livkid$v.plot, layout = jlayout)  

dev.off()


# order by tissue.order
dat.long$tissue <- factor(dat.long$tissue, levels = tissue.order)
pdf(file.path(outdir, "tissue_modules.examples.pdf"))
# BMAL1
print(PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Arntl" & experiment == "rnaseq"))) 
print(PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Nr1d1" & experiment == "rnaseq")))
print(PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Myh4" & experiment == "rnaseq")))
print(PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Pvalb" & experiment == "rnaseq")))
print(PlotGeneAcrossTissuesRnaseq(subset(dat.long, gene == "Tnnt3" & experiment == "rnaseq")))

dev.off()


# Cluster tissues to overlay later

pdf(file.path(outdir, "tissue_cluster.pdf"))
tfs <- GetTFs()
dat.mean <- subset(dat.long, gene %in% tfs & experiment == "rnaseq" & !tissue %in% filt.tiss) %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs))
mat <- dcast(dat.mean, formula = gene ~ tissue, value.var = "exprs.mean") 
rownames(mat) <- mat$gene; mat$gene <- NULL
hc <- hclust(dist(t(mat)))
plot(hc)
dev.off()

# top pairs
fits.count <- fits.best %>%
  group_by(model) %>%
  summarise(count = length(gene))
fits.count$n.rhyth <- sapply(fits.count$model, GetNrhythFromModel)
fits.count <- fits.count[order(fits.count$count, decreasing = TRUE), ]
fits.count$model <- factor(fits.count$model, level = unique(fits.count$model))
# plot top genes for n.rhyth >= 2
pdf(file.path(outdir, "top_multi_models.pdf"))
for (jmodel in fits.count$model[2:50]){
  if (jmodel != ""){
    PlotGenesInModel(subset(fits.best, model == jmodel), dat.complex, filt.tiss)
  }
}
dev.off()

# Alternative promoters ---------------------------------------------------


