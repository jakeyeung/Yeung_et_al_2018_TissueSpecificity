# Jake Yeung
# 2015-09-11
# total_variance.R
# Plot total variance of R.


# Functions ---------------------------------------------------------------

library(dplyr)
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FourierFunctions.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LoadArray.R")
source("scripts/functions/VarianceFunctions.R")

# Load data ---------------------------------------------------------------

# dat.long <- LoadArrayRnaSeq(fix.rik.xgene = TRUE)
# dat.array.noadj <- LoadArray(form = "long")
# dat.array.noadj$exprs <- dat.array.noadj$signal; dat.array.noadj$signal <- NULL
# dat.array.noadj$experiment <- "array.not_adj"
load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/dat.array.noadj.Robj", verbose=T)

load(file = "Robjs/dat.fit.Robj")

ref.gene <- "Nr1d1"

# dat.fit.relamp <- GetRelamp(dat.fit, max.pval = 1e-3)
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)

dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))

# Transform Fourier -------------------------------------------------------

# x <- subset(dat.array, gene == "Elovl3" & tissue == "Liver")$exprs
# x.per <- CalculatePeriodogram(x)
# PlotPeriodogram(x.per$freq, x.per$p.scaled)  # max is 4 hours

genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)

dat.array <- subset(dat.long, experiment == "array")

dat.sub <- subset(dat.array, gene %in% genes.exprs)

periods <- rep(48, 12) / seq(1, 12)  # 48/1, 48/2 ... 48/12

library(parallel)
dat.complexes <- lapply(periods, function(period, jdat, add.entropy.method){
  dat.comp <- TemporalToFrequencyDatLong(dat = jdat, period = period, add.entropy.method = add.entropy.method)
  dat.comp$period <- period
  return(dat.comp)
# }, jdat = dat.long, add.entropy.method = FALSE)
}, jdat = dat.sub, add.entropy.method = FALSE)

# save(dat.complexes, file = "Robjs/dat.complex.all_periods.array_only.Robj")

dat.complex.all_T <- do.call(rbind, dat.complexes)


# Explore noisy genes -----------------------------------------------------

jdat.4 <- GetTopGenesByPeriod(dat.complex.all_T, 4)
head(jdat.4)  # Saa2 in Liver

PlotGeneAcrossTissues(subset(dat.long, gene == "Saa2"))
PlotGeneAcrossTissues(subset(dat.array.noadj, gene == "Saa2"))


# Variance of each s  -----------------------------------------------------

dat.var.s <- dat.complex.all_T %>%
  group_by(period, tissue) %>%
  summarise(sum_sqr_mod = sum(Mod(exprs.transformed) ^ 2)) %>%
  mutate(period.factor = signif(period, digits = 3))

dat.var.s$period.factor <- factor(dat.var.s$period.factor, 
                           levels = sort(unique(dat.var.s$period.factor), decreasing = TRUE))


# Plot distribution of variances for each tissue --------------------------

ggplot(dat.var.s, aes(x = period.factor, y = sum_sqr_mod)) +  geom_bar(stat = "identity") + facet_wrap(~tissue) + 
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) + xlab("Period (hour)") + ylab("Sum of modulus squared across genes")


# Sum of all variances across all tissues ---------------------------------

dat.var.all <- dat.var.s %>%
  group_by(tissue) %>%
  summarise(sum_sqr_mod.total = sum(sum_sqr_mod))
dat.var.all$tissue <- factor(dat.var.all$tissue,
                             levels = dat.var.all$tissue[order(dat.var.all$sum_sqr_mod.total, decreasing = TRUE)])
# for ordering the facet_wrap plot across tissues by total variance
dat.var.s$tissue <- factor(dat.var.s$tissue,
                           levels = dat.var.all$tissue[order(dat.var.all$sum_sqr_mod.total, decreasing = TRUE)])

ggplot(dat.var.all, aes(x = tissue, y = sum_sqr_mod.total)) +  geom_bar(stat = "identity")  + 
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) + xlab("") + ylab("Sum of modulus squared across genes and s")

ggplot(dat.var.s, aes(x = tissue, y = sum_sqr_mod, fill = period.factor)) +  geom_bar(stat = "identity") + 
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) + xlab("") + ylab("Sum of modulus squared across genes and s")


# 24-h component normalized by total variance -----------------------------

dat.var.s1_adj <- dat.var.s %>%
  group_by(tissue) %>%
  mutate(s1_normalized = sum_sqr_mod / sum(sum_sqr_mod)) %>%
  filter(period == 24)
dat.var.s1_adj$tissue <- factor(dat.var.s1_adj$tissue,
                             levels = dat.var.s1_adj$tissue[order(dat.var.s1_adj$s1_normalized, decreasing = TRUE)])

# plot normalized spectral power
ggplot(dat.var.s1_adj, aes(x = tissue, y = s1_normalized)) + geom_bar(stat = "identity") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) + xlab("") + ylab("Normalized 24h spectral power")

# plot 24h spectral power
ggplot(OrderDecreasing(subset(dat.var.s, period == 24), "tissue", "sum_sqr_mod"), aes(x = tissue, y = sum_sqr_mod)) + geom_bar(stat = "identity") +
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) + xlab("") + ylab("24h spectral power (not normalized)")


# Why is amplitude between between Liver and Hypo so large in Nr1d --------

load("Robjs/dat.complex.maxexprs4.Robj", verbose = T)
head(dat.complex)
subset(dat.complex, gene == "Nr1d1")


# Plot distribution of genes yo -------------------------------------------

dat.complex <- dat.complex %>%
  group_by(tissue) %>%
  mutate(i.exprs.transformed = GetIndex(Mod(exprs.transformed)),
         i.exprs.adj = GetIndex(Mod(exprs.adj)))

ggplot(subset(dat.complex, i.exprs.transformed < 500), aes(x = i.exprs.transformed, y = Mod(exprs.transformed) ^ 2)) + 
  geom_line() + 
  facet_wrap(~tissue) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + xlab("Gene Index")

ggplot(subset(dat.complex, i.exprs.adj < 500), aes(x = i.exprs.adj, y = Mod(exprs.adj) ^ 2)) + 
  geom_line() + 
  facet_wrap(~tissue) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + xlab("Gene Index")


# dat.sub <- subset(dat.complex, tissue == "Hypo")
# dat.sub$i <- GetIndex(Mod(dat.sub$exprs.transformed))
# 
# # dat.sub$gene <- factor(as.character(dat.sub$gene), levels = dat.sub$gene[order(Mod(dat.sub$exprs.transformed), decreasing = TRUE)])
# ggplot(dat.sub[order(dat.sub$i), ], aes(x = i, y = Mod(exprs.transformed))) + 
#   geom_point() + 
#   facet_wrap(~tissue) + 
#   theme(axis.ticks = element_blank(), axis.text.x = element_blank())
