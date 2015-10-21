# 2015-10-20 
# Figures that will be part of paper to analyze this hogenesch dataset

remove.wfat <- TRUE

setwd("/home/yeung/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots_paper"

library(ggplot2)
library(hash)
library(dplyr)
library(reshape2)
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
ends <- c(12, 99, 288)
pdf(file.path(outdir, "PCA_scree.pdf"))
for (i in seq(length(starts))){
  print(PlotComponents2(pca.p.sum.long, starts[[i]], ends[[i]]))  # PCA 1,2,13,14,15,16, .. .20 are interesting
}
dev.off()

# plot first 50 modules across tissues and time
pdf(file.path(outdir, "PCA_modules.pdf"))
pcs <- seq(50); pcs <- paste("PC", pcs, sep = "")
for (jpc in pcs){
  print(ggplot(subset(pca.long, pc == jpc), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jpc))
}
dev.off()


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
# for ordering the facet_wrap plot across tissues by total variance
dat.var.s$tissue <- factor(dat.var.s$tissue,
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

pdf(file.path(outdir, "fourier_variance.pdf"))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(subset(dat.var.s), aes(x = tissue, y = sum_sqr_mod, fill = period.factor.cond)) +  geom_bar(stat = "identity") + 
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") + xlab("") + ylab("Spectral power") +
  scale_fill_manual(name = "Fourier component", drop=TRUE, values=cbPalette)
dev.off()

dat.var.s1_adj <- dat.var.s %>%
  group_by(tissue) %>%
  mutate(s1_normalized = sum_sqr_mod / sum(sum_sqr_mod))

dat.var.s1_adj.24 <- subset(dat.var.s1_adj, period.factor.cond == "24")
dat.var.s1_adj.24$tissue <- factor(dat.var.s1_adj.24$tissue,
                                levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])

pdf(file.path(outdir, "circadian_and_ultraidian_variance.pdf"))
# plot normalized spectral power
ggplot(dat.var.s1_adj.24, aes(x = tissue, y = s1_normalized)) + geom_bar(stat = "identity") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") + xlab("") + ylab("Normalized 24h spectral power") + 
  facet_wrap(~period.factor.cond)

dat.var.s1_adj.12 <- subset(dat.var.s1_adj, period.factor.cond == "12")
dat.var.s1_adj.12$tissue <- factor(dat.var.s1_adj.12$tissue,
                                   levels = dat.var.s1_adj.12$tissue[order(dat.var.s1_adj.12$s1_normalized, decreasing = TRUE)])

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
ggplot(dat.var.s1_adj, aes(x = period.factor, y = s1_normalized)) +  geom_bar(stat = "identity") + facet_wrap(~tissue) + 
  theme_bw(24) + theme(axis.text.x=element_text(size=11, angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") +
  facet_wrap(~tissue) +
  xlab("Fourier component (h)") + ylab("Normalized spectral power")
dev.off()

# Figure 1C ---------------------------------------------------------------
# Core clock genes

clockgenes <- GetClockGenes()

load(file = "Robjs/dat.fit.scan_periods.Robj")
load(file = "Robjs/dat.long.fixed_rik_genes.Robj")
dat.fit.periods <- subset(dat.fit.periods, gene %in% clockgenes)
if (remove.wfat){
  dat.fit.periods <- subset(dat.fit.periods, tissue != "WFAT")
}

dat.fit.periods$chi.sqr <- dat.fit.periods$ssq.residuals / dat.fit.periods$variance

dat.fit.periods <- subset(dat.fit.periods, gene != "Asb12")
ggplot(subset(dat.fit.periods), aes(x = period, y = chi.sqr, colour = gene, group = gene)) + geom_line() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")
# ggplot(subset(dat.fit.periods, ! gene %in% c("Asb12")), aes(x = period, y = ssq.residuals, colour = gene, group = gene)) + geom_line() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")


# Get minimum per gene ----------------------------------------------------

dat.fit.periods.min <- dat.fit.periods %>%
  group_by(gene, tissue) %>%
  do(GetMinPeriodSsqResiduals(.))

pdf(file.path(outdir, "fourier_across_periods.pdf"))
ggplot(subset(dat.fit.periods.min), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.min$period))/55) + geom_vline(xintercept=24, linetype="dotted") + scale_x_continuous(breaks=seq(20, 30, 1)) + 
  xlab("Period with minimum RSS") + ylab("Count") +
  theme_bw(24) + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       axis.text.x=element_text(size=11))
ggplot(subset(dat.fit.periods.min), aes(x = period)) + geom_histogram(binwidth = diff(range(dat.fit.periods.min$period))/55) + geom_vline(xintercept=24, linetype="dotted") + scale_x_continuous(breaks=seq(20, 30, 3)) + 
  facet_wrap(~tissue) +
  xlab("Period with minimum RSS") + ylab("Count") +
  theme_bw(24) + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       axis.text.x=element_text(size=11))
ggplot(dat.fit.periods.min, aes(x = period, y = chi.sqr, colour = gene)) + geom_point() + facet_wrap(~tissue) + geom_vline(xintercept=24, linetype="dotted")

# Observe Cry1 Aorta why it is period of 22
tiss <- "Aorta"; gen <- "Cry1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.min, tissue == tiss & gene == gen)$period
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
period.min <- subset(dat.fit.periods.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Mus"; gen <- "Nr1d1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Hypo"; gen <- "Nfil3"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Hypo"; gen <- "Wee1"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

tiss <- "Adr"; gen <- "Arntl"; exper="array"
dat.sub <- subset(dat.long, gene == gen & tissue == tiss)
period.min <- subset(dat.fit.periods.min, tissue == tiss & gene == gen)$period
PlotFitTwoPeriods(dat.sub, period1 = 24, period2 = period.min, tiss, gen, exper)

dev.off()
