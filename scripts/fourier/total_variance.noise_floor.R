# 2016-04-14
# Jake Yeung
# assess noise floor in total variance

rm(list=ls())

library(ggplot2)
library(dplyr)

remove.wfat <- TRUE

# Functions ---------------------------------------------------------------

source("scripts/functions/FourierFunctions.R")


# Load --------------------------------------------------------------------



load("Robjs/dat.var.filt.fixed.Robj", v=T)





# MAIN --------------------------------------------------------------------

# order factors
dat.var.filt.bytiss <- dat.var.filt %>%
  group_by(tissue) %>%
  summarise(var.temp = sum(var.temp), var.gene = sum(var.gene))
dat.var.filt$tissue <- factor(dat.var.filt$tissue, levels = dat.var.filt.bytiss$tissue[order(dat.var.filt.bytiss$var.temp, decreasing = T)])

dat.var.filt.sort <- dat.var.filt %>%
  group_by(tissue) %>%
  arrange(., desc(var.temp)) %>%
  mutate(var.temp.i = seq(length(var.temp)),
         var.temp.cum = cumsum(var.temp),
         var.temp.cum.norm = cumsum(var.temp) / (sum(var.temp) + sum(var.gene)))

ggplot(subset(dat.var.filt.sort, tissue != "WFAT"), aes(y = var.temp.cum.norm, x = var.temp.i)) + geom_line() + facet_wrap(~tissue) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("Cumulant of temporal variance (normalized by total variance)") + xlab("Gene index") + ggtitle("Contribution of temporal variance for each gene")


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
spectral.power <- ggplot(subset(dat.var.s, period.factor.cond != "Other"), aes(x = tissue, y = sum_sqr_mod, fill = period.factor.cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") + xlab("") + ylab("Spectral power") +
  scale_fill_manual(name = "Fourier component", drop=TRUE, values=cbPalette)
# add noise floor
periods <- sort(unique(dat.var.s$period))
noise.components <- periods[which(24 %% periods != 0)]

dat.var.s.noisefloor <- subset(dat.var.s, period %in% noise.components) %>%
  group_by(period.factor) %>%
  summarise(sum_sqr_mod.mean = mean(sum_sqr_mod), 
            sum_sqr_mod.var = var(sum_sqr_mod), 
            sum_sqr_mod.min = min(sum_sqr_mod), 
            sum_sqr_mod.max = max(sum_sqr_mod))

noise.floor <- mean(dat.var.s.noisefloor$sum_sqr_mod.mean)
spectral.power <- spectral.power + geom_hline(aes(yintercept = noise.floor), linetype="dotted")
print(spectral.power)

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

dat.var.s1_adj$tissue <- factor(dat.var.s1_adj$tissue,
                                levels = dat.var.s1_adj.24$tissue[order(dat.var.s1_adj.24$s1_normalized, decreasing = TRUE)])

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

m1.unadj <- ggplot(dat.var.s1_adj, aes(x = period.factor, y = sum_sqr_mod)) +  geom_bar(stat = "identity") + facet_wrap(~tissue) + 
  theme_bw(24) + theme(axis.text.x=element_text(size=11, angle=90,vjust = 0, hjust = 1), 
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "black"),
                       legend.position="bottom") +
  facet_wrap(~tissue) +
  xlab("Fourier component (h)") + ylab("Spectral power")
print(m1.unadj)


# N genes by amp ----------------------------------------------------------


# plot log amplitude and number of genes above amplitude

load("Robjs/dat.fit.Robj", v=T); dat.fit.24 <- dat.fit
load("Robjs/dat.fit.period.12.Robj", v=T); dat.fit.12 <- dat.fit
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)

source("scripts/functions/PlotGeneAcrossTissues.R", v=T)

if (remove.wfat){
  dat.long <- subset(dat.long, tissue != "WFAT")
  dat.fit.24 <- subset(dat.fit.24, tissue != "WFAT")
  dat.fit.12 <- subset(dat.fit.12, tissue != "WFAT")
}

PlotGeneAcrossTissues(subset(dat.long, gene == "Nr1d1"))

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

ggplot(dat.fit.24.ngenes.thres, aes(x = 2 * amp.thres, y = n.genes)) + geom_line() + 
  facet_wrap(~tissue) + theme_bw(24) + ggtitle(paste("24 hour rhythms. Pval >", pval.cutoff)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Amplitude (fold change peak/trough)") + ylab("Number of Genes")
ggplot(dat.fit.12.ngenes.thres, aes(x = 2 * amp.thres, y = n.genes)) + geom_line() + 
  facet_wrap(~tissue) + theme_bw(24) + ggtitle(paste("12 hour rhythms. Pval >", pval.cutoff)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Amplitude (fold change peak/trough)") + ylab("Number of Genes")

# combine the two
dat.fit.24.ngenes.thres$rhyth <- as.factor(24)
dat.fit.12.ngenes.thres$rhyth <- as.factor(12)
dat.fit.ngenes.thres <- rbind(dat.fit.24.ngenes.thres, dat.fit.12.ngenes.thres)

ggplot(dat.fit.ngenes.thres, aes(x = 2 * amp.thres, y = n.genes + 1, colour = rhyth)) + geom_line() + 
  facet_wrap(~tissue) + theme_bw(24) + ggtitle(paste("12 and 24 hour rhythms. Pval >", pval.cutoff)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Amplitude (fold change peak/trough)") + ylab("Number of Genes") + xlim(c(0, 5)) + 
  scale_y_log10()

dat.test <- subset(dat.fit.24, pval < 0.01 & tissue == "Liver")

n.genes.amp.thres <- NGenesByAmp(dat.test$amp, amp.thres)

# order by total genes
ngenes.sum <- dat.fit.ngenes.thres %>%
  group_by(tissue) %>%
  summarise(total = sum(n.genes)) %>%
  arrange(desc(total))
dat.fit.ngenes.thres$tissue <- factor(as.character(dat.fit.ngenes.thres$tissue), levels = ngenes.sum$tissue)
ggplot(subset(dat.fit.ngenes.thres, rhyth == 24), aes(x = 2 * amp.thres, y = n.genes, colour = tissue)) + geom_line() + 
  theme_bw(24) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Log2 Fold Change") + ylab("# Genes") + xlim(c(0, 5)) + 
  scale_y_log10() + 
  geom_vline(xintercept = 2, linetype = "dotted")

# How much variance is explained by our modules ---------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
source("scripts/functions/PlotFunctions.R")
genes.tw <- as.character(subset(fits.best, n.rhyth >= 8)$gene)
genes.adr <- as.character(subset(fits.best, model == "Adr")$gene)
genes.liver <- as.character(subset(fits.best, model == "Liver")$gene)
genes.bfat <- as.character(subset(fits.best, model == "BFAT")$gene)
genes.mus <- as.character(subset(fits.best, model == "Mus")$gene)

genes.lst <- list("TissueWide"=genes.tw, "Adr"=genes.adr, "Liver"=genes.liver, "BFAT"=genes.bfat, "Mus"=genes.mus)
frac.var.lst <- list()

for (gene.group in names(genes.lst)){
  print(gene.group)
  genes.sub <- genes.lst[[gene.group]]
  frac.var.lst[[gene.group]] <- FracVarByGeneList(dat.var.s, dat.var.filt, dat.complex.all_T, 
                                                  genes.lst[[gene.group]], period.for.sort = 24, jlab = gene.group)
}
frac.var.lst <- do.call(rbind, frac.var.lst)

ggplot(frac.var.lst, aes(x = tissue, y = var.norm, fill = as.factor(period))) + geom_bar(stat="identity", position=position_dodge()) + facet_wrap(~label)

# SVD on residual genes ---------------------------------------------------

load("Robjs/dat.complex.fixed_rik_genes.Robj")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")

if (remove.wfat){
  dat.complex <- subset(dat.complex, tissue != "WFAT")
}

fits.res <- subset(fits.best, n.rhyth < 7 & n.rhyth > 2)
fits.res <- subset(fits.best, n.rhyth == 2)
genes.res <- as.character(fits.res$gene)

s.res <- SvdOnComplex(subset(dat.complex, gene %in% genes.res), value.var = "exprs.transformed")
eigens.res <- GetEigens(s.res, period = 24, comp = 3, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.res$u.plot, eigens.res$v.plot, layout = jlayout)  



