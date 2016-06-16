# 2016-06-15
# Jake Yeung

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

plotdir <- "plots/liver_kidney_math710"
dir.create(plotdir)

library(dplyr)
library(ggplot2)
setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FitRhythmicAcrossPeriods.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/ListFunctions.R")

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
# Load --------------------------------------------------------------------

dat <- LoadLivKid()
foutpath <- "Robjs/liver_kidney/dat.fit.periods.genome_wide.Robj"
load(foutpath, verbose = T)

indir <- "/home/yeung/projects/tissue-specificity/Robjs/bayes_factors"
inf <- list.files(indir)
fits.long <- expandingList()
for (f in inf){
  print(f)
  method <- strsplit(f, split = "\\.")[[1]][[5]]
  fpath <- file.path(indir, f)
  load(fpath)
  fits.all.long$method <- method
  fits.long$add(fits.all.long)
}
fits.long <- fits.long$as.list()
fits.long <- do.call(rbind, fits.long)

# Process lowly expressed genes -------------------------------------------

# x <- dat$exprs[which(dat$exprs > min(dat$exprs))]
# x <- sample(x, 0.01 * length(x))
# plot(density(x)); abline(v=1)
# source("scripts/functions/MixtureModelFunctions.R")
# cutoff <- FindCutoff(x, lambdas = c(0.01, 0.99), mus = c(-1, 4), k = 2, show.fig = TRUE)$maximum
cutoff <- -0.187448

dat.mean <- dat %>%
  group_by(gene) %>%
  summarise(exprs.mean = mean(exprs))

dat.mean.bytiss <- dat %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs))

genes.filt <- as.character(subset(dat.mean, exprs.mean < cutoff)$gene)

# FILTER MAX AND MIN PERIODS 
# dat.fit.periods.genome_wide <- subset(dat.fit.periods.genome_wide)

# FILTER LOWLY EXPRSED GENES
dat <- subset(dat, ! gene %in% genes.filt)
fits.long <- subset(fits.long, ! gene %in% genes.filt)
# get best model
fits.long.filt <- fits.long %>%
  group_by(gene, method) %>%
  filter(weight.raw == min(weight.raw))
# frequency space
omega <- 2 * pi / 24
dat.freq <- dat %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))

# Param sweep -------------------------------------------------------------

dat.fit.periods.genome_wide.min <- dat.fit.periods.genome_wide %>%
  group_by(gene, tissue) %>%
  filter(period == period[which.min(ssq.residuals)])

xscale_periods <- seq(12, 28, 4)
plot.periods.all <- ggplot(dat.fit.periods.genome_wide.min, aes(x = period)) + 
  geom_histogram(binwidth = 0.5) + 
  geom_vline(xintercept=c(12, 24), linetype="dotted") + 
  scale_x_continuous(breaks=xscale_periods) + 
  xlab("Period with best fit [h]") + ylab("Number of genes") +
  theme_bw() + 
  theme(aspect.ratio=1,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.text.x = element_text(size=15),
          axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24))

# Variance of each s  -----------------------------------------------------

load("Robjs/liver_kidney/dat.complex.all_T.livkid.Robj")

# multiply by 2 to get amplitude
dat.complex.all_T <- dat.complex.all_T %>%
  group_by(gene, tissue) %>%
  mutate(exprs.transformed = exprs.transformed * 2)

dat.var.s <- dat.complex.all_T %>%
  group_by(period, tissue) %>%
  summarise(sum_sqr_mod = sum(Mod(exprs.transformed) ^ 2)) %>%
  mutate(period.factor = signif(period, digits = 3))

dat.var.s$period.factor <- factor(dat.var.s$period.factor, 
                                  levels = sort(unique(dat.var.s$period.factor), decreasing = TRUE))
dat.var.s <- dat.var.s %>%
  group_by(tissue) %>%
  mutate(sum_sqr_mod.norm = sum_sqr_mod / sum(sum_sqr_mod))

plot.fourier <- ggplot(dat.var.s, aes(x = period.factor, y = sum_sqr_mod)) +  geom_bar(stat = "identity") + facet_wrap(~tissue) + 
  theme_bw(24) + theme(aspect.ratio=1, 
                       axis.text.x=element_text(angle=90,vjust = 0, hjust = 1),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank()) + 
  xlab("Period [h]") + ylab(expression(paste("Spectral power [", log[2], RPKM^2, "]", sep="")))
plot.fourier.norm <- ggplot(dat.var.s, aes(x = period.factor, y = sum_sqr_mod.norm)) +  geom_bar(stat = "identity") + facet_wrap(~tissue) + 
  theme_bw(24) + theme(aspect.ratio=1, 
                       axis.text.x=element_text(angle=90,vjust = 0, hjust = 1),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank()) + 
  xlab("Period [h]") + ylab("Normalized spectral power")


# Nconds Summary ----------------------------------------------------------

fits.sum <- fits.long.filt %>%
  group_by(method, model) %>%
  summarise(count = length(model))

fits.sum$n.params <- sapply(as.character(fits.sum$model), function(m){
  n.tiss <- 2
  s <- strsplit(m, ";")[[1]]
  n.params <- n.tiss + 2 * length(s)
  return(paste0("(p=", n.params, ")"))
})

fits.sum$model.label <- mapply(function(model, n.param){
  if (model == "") model <- "Flat"
  paste(model, n.param, sep = "\n")
}, fits.sum$model, fits.sum$n.params)

fits.sum$model.label <- factor(fits.sum$model.label, levels = fits.sum$model.label[1:5])

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#D55E00", "#CC79A7", "#F0E442", "#009E73", "#0072B2")
plot.counts <- ggplot(fits.sum, aes(x = model.label, fill = method, y = count)) + geom_bar(stat = "identity", position = "dodge", width=0.5) + theme_bw(24) + 
  scale_fill_manual(values=cbPalette, labels = c("AIC", "BIC", "g=1000", "HyperG", "Zellner-Siow")) + xlab("") + 
  theme(aspect.ratio = 1, axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())

plot.counts.noeb <- ggplot(subset(fits.sum, method != "eb"), aes(x = model.label, fill = method, y = count)) + geom_bar(stat = "identity", position = "dodge", width=0.5) + theme_bw(24) + 
  scale_fill_manual(values=cbPalette, labels = c("AIC", "BIC", "HyperG", "Zellner-Siow")) + xlab("") + 
  theme(aspect.ratio = 1, axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())


# Discrepancies -----------------------------------------------------------

fits.disc <- fits.long.filt %>%
  group_by(gene) %>%
  filter(length(unique(model)) > 1)
fits.disc$amp <- sapply(fits.disc$param.list, function(p) GetAvgAmpFromParams(params = p, by.model = FALSE))
fits.disc$phase <- sapply(fits.disc$param.list, function(p) GetPhaseFromParams(params = p, by.model = FALSE))
fits.disc$phase.diff <- sapply(fits.disc$param.list, function(p) GetMaxPhaseDiffFromParams(params = p, by.model=FALSE))

fits.disc <- fits.disc[order(fits.disc$amp, decreasing=TRUE), ]

# SVD Downstrream ---------------------------------------------------------

# also print everything

pdf(file.path(plotdir, "math710plots.pdf"))
print(plot.periods.all)
print(plot.periods.all + facet_wrap(~tissue) + theme(strip.text.x = element_text(size = 12)))
print(plot.fourier)
print(plot.fourier.norm)
print(plot.counts)
dev.off()

pdf(file.path(plotdir, "liverkidneysvd.pdf"))
# Compare Kidney,Liver model
for (jmeth in c("AIC", "hyperg", "BIC", "zf", "eb")){
  for (jtiss in list("Kidney,Liver")){
    genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% jtiss)$gene)
    for (i in c(1)){
      s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
      eigens1 <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
      eigens2 <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE, jtitle = jmeth)
      jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
      multiplot(eigens1$u.plot, eigens1$v.plot, layout = jlayout)
      multiplot(eigens2$u.plot, eigens2$v.plot, layout = jlayout)
    }
  }  
}
dev.off()

pdf(file.path(plotdir, "liverkidneysepparams.pdf"))
# Compare Kidney,Liver model
for (jmeth in c("AIC", "hyperg", "BIC", "zf", "eb")){
  for (jtiss in list("Kidney;Liver")){
    genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% jtiss)$gene)
    for (i in c(1, 2)){
      s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
      eigens1 <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
      eigens2 <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE, jtitle = jmeth)
      jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
      multiplot(eigens1$u.plot, eigens1$v.plot, layout = jlayout)
      multiplot(eigens2$u.plot, eigens2$v.plot, layout = jlayout)
    }
  }  
}
dev.off()

pdf(file.path(plotdir, "genomewide.pdf"))
# Compare Kidney,Liver model
      s <- SvdOnComplex(subset(dat.freq), value.var = "exprs.transformed")
      eigens1 <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
      eigens2 <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE, jtitle = jmeth)
      jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
      multiplot(eigens1$u.plot, eigens1$v.plot, layout = jlayout)
      multiplot(eigens2$u.plot, eigens2$v.plot, layout = jlayout)
dev.off()

for (jtiss in list("Liver", "Kidney")){
  outf <- paste0(jtiss, ".pdf")
  pdf(file.path(plotdir, outf))
  # Compare Kidney,Liver model
  for (jmeth in c("AIC", "hyperg", "BIC", "zf", "eb")){
    genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% jtiss)$gene)
    for (i in c(1)){
      s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
      eigens1 <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
      eigens2 <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE, jtitle = jmeth)
      jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
      multiplot(eigens1$u.plot, eigens1$v.plot, layout = jlayout)
      multiplot(eigens2$u.plot, eigens2$v.plot, layout = jlayout)
    }
  }
  dev.off()
}


# plot examples
jgenes <- c("Arntl", "Nr1d1", "Npas2", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Gm15459", "Hspa8", "Dbp", "Slc44a1", "Osgin1", "Gm11128")
for (jgene in jgenes){
  print(PlotGeneAcrossTissues(subset(dat, gene == jgene), make.pretty = TRUE))
  print(PlotGeneAcrossTissues(subset(dat, gene == jgene), make.pretty = TRUE, do.facet.wrap = FALSE))
  print(ggplot(subset(fits.long, gene == jgene), aes(x = method, y = weight, fill = model)) + geom_bar(stat = "identity", position = "dodge") + theme_bw() + ggtitle(jgene) + 
    scale_x_discrete(labels = c("AIC", "BIC", "g=1000", "HyperG", "Zellner-Siow")) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())) + 
    xlab("") + ylab("Probability")
  print(PlotPeriodogramLong(subset(dat, gene == jgene)))
}
dev.off()

# Plot outputs ------------------------------------------------------------

# 12 hours

jsub <- subset(dat.fit.periods.genome_wide.min, abs(period - 12) <= 2 & int > 2.5 & amp > 0.25)
jsub <- jsub[order(jsub$amp, decreasing = TRUE), ]
jsub <- jsub[order(jsub$pval, decreasing = FALSE), ]
jgenes <- as.character(jsub$gene)

pdf(file.path(plotdir, "ultradiangenes.pdf"))
for (jgene in jgenes){
  print(PlotGeneAcrossTissues(subset(dat, gene == jgene), make.pretty = TRUE))
  print(PlotPeriodogramLong(subset(dat, gene == jgene)))
}
dev.off()

# 
# jgene <- "Gm15459"
# jgene <-"Hspa8"
# jgene <-"Dbp"
# PlotGeneAcrossTissues(subset(dat, gene == jgene), make.pretty = TRUE)
# PlotPeriodogramLong(subset(dat, gene == jgene))

# PlotGeneAcrossTissues((subset(dat, gene == "Gm15459")))
# 
# PlotGeneAcrossTissues((subset(dat, gene == "Scd1")))
# PlotGeneAcrossTissues((subset(dat, gene == "Lars2")))
# PlotGeneAcrossTissues((subset(dat, gene == "G6pc")))
# PlotGeneAcrossTissues((subset(dat, gene == "Hspa13")))
# PlotGeneAcrossTissues((subset(dat, gene == "Fnbp4")))

