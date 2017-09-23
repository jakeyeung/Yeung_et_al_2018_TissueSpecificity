# Jake Yeung
# Date of Creation: 2017-09-22
# File: ~/projects/tissue-specificity/scripts/primetime_figures/tissue_specificity_paper_rerun_PCA.R
# Rerun PCA

remove.wfat <- TRUE

tissue.order <- c('Liver','BFAT','Kidney','Lung','Adr','Mus','Heart','Aorta','Hypo','Cere','BS')

setwd("/home/yeung/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots_paper_rerun_PCA"
dir.create(outdir)

library(ggplot2)
library(hash)
library(dplyr)
library(reshape2)
library(gplots)
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
source("scripts/functions/LoadActivitiesLong.R")


# Figure 1A Present high level variation BETWEEN TISSUES and WITHIN TISSUES (time) ------------------
# scripts from pca_adjusted_microarray.label_variance.R

# Load data, log transform

log2.transform <- TRUE
# load("Robjs/dat.array.adj.primetime.Robj", verbose=T)  # no -Infs

dat <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt",
                           remove.negs = TRUE, fix.rik.xgene = TRUE)
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

plot.i <- 1
pdf(file.path(outdir, paste0(plot.i, ".component_vs_component.pdf")), useDingbats=FALSE)
plot.i <- plot.i + 1
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


pdf(file.path(outdir, paste0(plot.i, ".PCA_scree.pdf")))
plot.i <- plot.i + 1

for (i in seq(length(starts))){
  print(PlotComponents2(pca.p.sum.long, starts[[i]], ends[[i]]))  # PCA 1,2,13,14,15,16, .. .20 are interesting
}
dev.off()

# plot first 50 modules across tissues and time

# order pca.p.sum.long tissues to same as Fourier components: hardcode
pca.long$tissue <- factor(pca.long$tissue, levels = tissue.order)

pdf(file.path(outdir, paste0(plot.i, ".PCA_modules.pdf")))
plot.i <- plot.i + 1

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
                 theme(aspect.ratio=1,
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
                 scale_x_continuous(breaks=c(18, 42, 64))
jpc <- "PC13"
pc2 <- ggplot(subset(pca.long, pc == jpc), aes(x = time, y = loading)) + geom_point() + geom_line() + facet_wrap(~tissue) + ggtitle(jpc) + 
  xlab("CT") +
  ylab("PC loading") + 
  theme_bw(24) + 
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
  scale_x_continuous(breaks=c(18, 42, 64))


multiplot(scree1, scree2, pc1, pc2, cols = 2)  


# Figure 1A var by tissue ---------------------------------------------------
# from pca_by_tissue.R
# Plot contribution of temporal variance for each gene across tissues
# load("Robjs/dat.var.filt.Robj")
load("Robjs/dat.var.filt.fixed.Robj", v=T)


# order factors
dat.var.filt.bytiss <- dat.var.filt %>%
  group_by(tissue) %>%
  summarise(var.temp = sum(var.temp), var.gene = sum(var.gene))
dat.var.filt$tissue <- factor(dat.var.filt$tissue, levels = dat.var.filt.bytiss$tissue[order(dat.var.filt.bytiss$var.temp, decreasing = T)])

pdf(file.path(outdir, paste0(plot.i, ".variance_by_tissue.pdf")))
plot.i <- plot.i + 1
dat.var.filt.sort <- dat.var.filt %>%
  group_by(tissue) %>%
  arrange(., desc(var.temp)) %>%
  mutate(var.temp.i = seq(length(var.temp)),
         var.temp.cum = cumsum(var.temp),
         var.temp.cum.norm = cumsum(var.temp) / (sum(var.temp) + sum(var.gene)))

ggplot(subset(dat.var.filt.sort, tissue != "WFAT"), aes(y = var.temp.cum.norm, x = var.temp.i)) + geom_line() + facet_wrap(~tissue) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("Cumulant of temporal variance (normalized by total variance)") + xlab("Gene index") + ggtitle("Contribution of temporal variance for each gene")

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

pdf(file.path(outdir, paste0(plot.i, ".fourier_variance.pdf")))
plot.i <- plot.i + 1
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
# dat.var.s$tissue <- factor(dat.var.s$tissue,
#                            levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])
dat.var.s$tissue <- factor(dat.var.s$tissue,
                           levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$sum_sqr_mod, decreasing = TRUE)])
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

pdf(file.path(outdir, paste0(plot.i, ".circadian_and_ultraidian_variance.pdf")))
plot.i <- plot.i + 1
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
pdf(file.path(outdir, paste0(plot.i, ".normalized_fourier_variance.pdf")))
plot.i <- plot.i + 1

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
