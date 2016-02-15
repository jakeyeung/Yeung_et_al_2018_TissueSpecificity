# 2015-12-01
# Split variance robustly so it is independent of number of tissues analyzed

library(ggplot2)
library(hash)
library(dplyr)
library(reshape2)

# Functions ---------------------------------------------------------------


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
source("scripts/functions/PlotFunctions.R")


GetTissueTemporalVar <- function(dat.gene, mean.global){
  # Given a gene over time, calculate its temporal variance (how the gene changes about its own mean) as well as
  # its gene mean variance (how the mean is different about its global mean)
  n.samp <- length(dat.genesub$exprs)
  dat.var.temp <- var(dat.gene$exprs)
  
  if (is.numeric(mean.global)){
    mean.all <- mean.global
  } else {
    # assume is hash
    tiss <- as.character(dat.gene$tissue[[1]])
    mean.all <- mean.global[[tiss]]
  }
  dat.var.gene <- sum((rep(mean(dat.gene$exprs), n.samp) - mean.all) ^ 2) / (n.samp - 1)
  return(data.frame(var.temp = dat.var.temp, var.gene = dat.var.gene))
}


# Load dat ----------------------------------------------------------------


dat <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt", 
                           remove.negs = TRUE, fix.rik.xgene = TRUE)

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)

dat.sub <- subset(dat.long, experiment == "array")


# Look at Liver -----------------------------------------------------------

globalmeans <- dat.sub %>%
  group_by(tissue) %>%
  summarise(mean.global = mean(exprs))

hash.means <- hash(as.character(globalmeans$tissue), globalmeans$mean.global)

dat.liv <- subset(dat.sub, tissue == "Liver")

# global mean
dat.liv.globalmean <- mean(dat.liv$exprs)

# variance of gene like DBP

jgene <- "Dbp"

dat.genesub <- subset(dat.liv, gene == jgene)

dat.var.total <- sum((dat.genesub$exprs - dat.liv.globalmean) ^ 2) / (length(dat.genesub$exprs) - 1)

# break it down into temporal variance and tissue variance
dat.var.temp <- sum((dat.genesub$exprs - mean(dat.genesub$exprs)) ^ 2) / (length(dat.genesub$exprs) - 1)
dat.var.temp <- var(dat.genesub$exprs)
n.samp <- nrow(dat.genesub)
dat.var.gene <- sum((rep(mean(dat.genesub$exprs), n.samp) - dat.liv.globalmean) ^ 2) / (n.samp - 1)


# Temporal versus mean variance genomewide --------------------------------

dat.var <- dat.sub %>%
  group_by(gene, tissue) %>%
  do(GetTissueTemporalVar(., hash.means))


# Plot summary ------------------------------------------------------------

dat.var.bytiss <- dat.var %>%
  group_by(tissue) %>%
  summarise(var.temp = sum(var.temp), var.gene = sum(var.gene))

dat.var.melt <- melt(dat.var.bytiss, variable.name = "variance.type", value.name = "variance", id.vars = c("tissue"))

dat.var.frac <- dat.var %>%
  group_by(tissue) %>%
  summarise(var.temp.frac = sum(var.temp) / (sum(var.temp) + sum(var.gene)))

dat.plot <- subset(dat.var.melt, tissue != "WFAT")
dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity") + facet_wrap(~variance.type)

dat.plot <- subset(dat.var.melt, tissue != "WFAT")
dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity")

dat.plot <- subset(dat.var.melt, variance.type == "var.gene" & tissue != "WFAT")
dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity") + facet_wrap(~variance.type)

dat.plot <- subset(dat.var.melt, variance.type == "var.temp" & tissue != "WFAT")
dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity") + facet_wrap(~variance.type)


# Why does Cere have large var.gene? --------------------------------------

ggplot(dat.long, aes(x = tissue, y = exprs)) + geom_boxplot() + facet_wrap(~experiment)

# order by total var.gene
dat.var$tissue <- factor(dat.var$tissue, levels = dat.var.bytiss$tissue[order(dat.var.bytiss$var.gene, decreasing = T)])
ggplot(dat.var, aes(x = var.gene)) + geom_histogram(binwidth = 5) + facet_wrap(~tissue)

# order by total var.temp

# How many genes are expressed above background? --------------------------

cutoff <- 4

dat.rnaseq <- subset(dat.long, experiment == "rnaseq")

dat.max <- dat.rnaseq %>%
  group_by(tissue, gene) %>%
  filter(max(exprs) > cutoff) %>%
  group_by(tissue) %>%
  summarise(n.exprs = length(exprs))

dat.min <- dat.rnaseq %>%
  group_by(tissue, gene) %>%
  filter(max(exprs) <= cutoff) %>%
  group_by(tissue) %>%
  summarise(n.exprs = length(exprs))

dat.max <- dat.rnaseq %>%
  group_by(tissue, gene) %>%
  summarise(exprs.max = max(exprs))

jkeys <- paste(as.character(dat.max$tissue), as.character(dat.max$gene), sep = ";")
jvals <- (dat.max$exprs.max >= cutoff)
hash.max <- hash(jkeys, jvals)

# Filter for genes expressed above background -----------------------------

dat.sub$is.exprs <- mapply(function(tiss, jgene) hash.max[[paste(tiss, jgene, sep = ";")]], as.character(dat.sub$tissue), as.character(dat.sub$gene))

# get global means across tissues after filtering for background
dat.means.filt <- subset(dat.sub, is.exprs == TRUE) %>%
  group_by(tissue) %>%
  summarise(exprs.mean = mean(exprs))

hash.means.filt <- hash(as.character(dat.means.filt$tissue), dat.means.filt$exprs.mean)

dat.var.filt <- subset(dat.sub, is.exprs == TRUE) %>%
  group_by(gene, tissue) %>%
  do(GetTissueTemporalVar(., mean.global = hash.means.filt))

# Plot summary ------------------------------------------------------------

dat.var.filt.bytiss <- dat.var.filt %>%
  group_by(tissue) %>%
  summarise(var.temp = sum(var.temp), var.gene = sum(var.gene))

dat.var.filt.melt <- melt(dat.var.filt.bytiss, variable.name = "variance.type", value.name = "variance", id.vars = c("tissue"))

dat.var.filt.frac <- dat.var.filt %>%
  group_by(tissue) %>%
  summarise(var.temp.frac = sum(var.temp) / (sum(var.temp) + sum(var.gene)))

dat.plot <- subset(dat.var.filt.melt, tissue != "WFAT")
dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity") + facet_wrap(~variance.type)

dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity")

dat.plot <- subset(dat.var.filt.melt, tissue != "WFAT")
dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity") + facet_wrap(~variance.type)

dat.plot <- subset(dat.var.filt.melt, variance.type == "var.gene" & tissue != "WFAT")
dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity") + facet_wrap(~variance.type)

dat.plot <- subset(dat.var.filt.melt, variance.type == "var.temp" & tissue != "WFAT")
dat.plot <- OrderDecreasing(dat.plot, jfactor = "tissue", jval = "variance")
ggplot(dat.plot, aes(x = tissue, y = variance, fill = variance.type)) + 
  geom_bar(stat="identity") + facet_wrap(~variance.type)



# order by total var.gene
dat.var.filt$tissue <- factor(dat.var.filt$tissue, levels = dat.var.filt.bytiss$tissue[order(dat.var.filt.bytiss$var.gene, decreasing = T)])
ggplot(subset(dat.var.filt, tissue != "WFAT"), aes(x = log10(var.gene))) + geom_histogram(binwidth = 0.1) + facet_wrap(~tissue)

dat.var.filt$tissue <- factor(dat.var.filt$tissue, levels = dat.var.filt.bytiss$tissue[order(dat.var.filt.bytiss$var.gene, decreasing = T)])
ggplot(subset(dat.var.filt, tissue != "WFAT"), aes(x = log10(var.gene))) + geom_histogram(binwidth = 0.4) + facet_wrap(~tissue)

dat.var.filt$tissue <- factor(dat.var.filt$tissue, levels = dat.var.filt.bytiss$tissue[order(dat.var.filt.bytiss$var.temp, decreasing = T)])
ggplot(subset(dat.var.filt, tissue != "WFAT"), aes(x = log10(var.temp))) + geom_histogram(binwidth = 0.05) + facet_wrap(~tissue)


# Split temporal variance into fourier components -------------------------

source("scripts/functions/SvdFunctions.R")
dat.isexprs <- subset(dat.sub, is.exprs == TRUE & tissue != "WFAT")

periods <- rep(48, 12) / seq(1, 12)  # 48/1, 48/2 ... 48/12

library(parallel)
dat.complexes <- lapply(periods, function(period, jdat, add.entropy.method){
  dat.comp <- TemporalToFrequencyDatLong(dat = jdat, period = period, add.entropy.method = add.entropy.method)
  dat.comp$period <- period
  return(dat.comp)
  # }, jdat = dat.long, add.entropy.method = FALSE)
}, jdat = dat.isexprs, add.entropy.method = FALSE)

dat.complexes <- do.call(rbind, dat.complexes)

dat.complexes$exprs.transformed <- as.complex(dat.complexes$exprs.transformed)

dat.complexes.filt <- subset(dat.complexes, !is.na(exprs.transformed))

dat.complexes.genomewide <- dat.complexes.filt %>%
  group_by(tissue) %>%
  summarise(sum.mod.sqr = sum(Mod(exprs.transformed) ^ 2))

# this should be equivalent to total temporal variance
dat.complexes.genomewide <- OrderDecreasing(dat.complexes.genomewide, jfactor = "tissue", jval = "sum.mod.sqr")
ggplot(dat.complexes.genomewide, aes(x = tissue, y = sum.mod.sqr)) + geom_bar(stat = "identity")

# plot variance by periods
dat.complexes.byperiod <- dat.complexes.filt %>%
  group_by(tissue, period) %>%
  summarise(sum.mod.sqr = sum(Mod(exprs.transformed) ^ 2)) %>%
  mutate(period.factor = signif(period, digits = 3))
dat.complexes.byperiod$period.factor <- factor(dat.complexes.byperiod$period.factor, 
                                  levels = sort(unique(dat.complexes.byperiod$period.factor), decreasing = TRUE))

dat.complexes.byperiod$tissue <- factor(as.character(dat.complexes.byperiod$tissue), 
                                        levels = dat.complexes.genomewide[order(dat.complexes.genomewide$sum.mod.sqr, decreasing = TRUE), ]$tissue)
ggplot(dat.complexes.byperiod, aes(x = period.factor, y = sum.mod.sqr)) + geom_bar(stat = "identity") + facet_wrap(~tissue)

# normalize to 1
dat.complexes.byperiod <- dat.complexes.byperiod %>%
  group_by(tissue) %>%
  mutate(sum.mod.sqr.norm = sum.mod.sqr / sum(sum.mod.sqr))

# order by decreasing 24 hour
# plot variance by periods normalized
dat.complexes.byperiod$tissue <- factor(as.character(dat.complexes.byperiod$tissue),
                                        levels = subset(dat.complexes.byperiod, 
                                                        period.factor == "24")[order(subset(dat.complexes.byperiod, 
                                                                                            period.factor == "24")$sum.mod.sqr.norm, 
                                                                                     decreasing = TRUE), ]$tissue)
ggplot(dat.complexes.byperiod, aes(x = period.factor, y = sum.mod.sqr.norm)) + geom_bar(stat = "identity") + facet_wrap(~tissue)
