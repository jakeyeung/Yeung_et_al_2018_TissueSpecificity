# Afte rrunning find_alternative_first_exons.R, I saved my objects so we can explore them later.


# Functions ---------------------------------------------------------------

setwd("~/projects/tissue-specificity/")

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/LoadAndHandleData.R")
# source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
# source("scripts/functions/MakeCluster.R")
# source("scripts/functions/ReadListToVector.R")
# source("scripts/functions/GrepRikGenes.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/TissueSpecificRhythmicsFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/LoadSitecounts.R")
source("scripts/functions/LoadLong.R")
source("scripts/functions/FitMotifAmp.R")

library(ggplot2)
library(dplyr)
library(mixtools)
library(gplots)
library(parallel)
# library(biglm)
library(reshape2)
library(hash)
# Functions ---------------------------------------------------------------



# Load my data ------------------------------------------------------------

# tissues <- c("Liver", "Mus", "Kidney", "Heart", "Lung")
tissues <- c("Liver", "Mus", "Kidney", "Heart", "Lung", "Cere")
encodetissues <- c("Liver", "SkeletalMuscle", "Kidney", "Heart", "Lung", "Cerebellum")
# encodetissues <- c("Liver", "SkeletalMuscle", "Kidney", "Heart", "Lung")
tissuesdic <- setNames(object = encodetissues, nm = tissues)

load(file = "results/alternative_exon_usage/cov.normreads.filt.rhyth.Robj")  # cov.normreads.filt.rhyth
load(file = "Robjs/fits.relamp.Robj")
ka.long <- LoadLong()

# ampdic.motif <- setNames(object = dat.rhyth.sub.motif$relamp, nm = paste(tissuesdic[as.character(dat.rhyth.sub.motif$tissue)], dat.rhyth.sub.motif$gene, sep = ';'))
vals <- fits.relamp$relamp
keys <- paste(tissuesdic[as.character(fits.relamp$tissue)], fits.relamp$gene, sep = ';')
ampdic.motif <- hash(keys, vals)

# Explore -----------------------------------------------------------------

head(cov.normreads.filt.rhyth)

# run model on full dataset
fit.afe <- cov.normreads.filt.rhyth %>%
  filter(!(is.na(rhythmic.or.not))) %>%
  group_by(transcript, gene) %>%
  do(FitRhythNonRhyth(jdf = .)) %>%
  filter(!is.na(pval))

# show top hits
(head(data.frame(fit.afe[order(fit.afe$pval), ]), n = 50))
# (head(data.frame(fit.afe[order(fit.afe$coef, decreasing = TRUE), ]), n = 50))

# summarize by choosing the top for each gene
fit.afe.summary <- fit.afe %>%
  group_by(gene) %>%
  do(SubsetMinPval(jdf = .))

(head(data.frame(fit.afe.summary[order(fit.afe.summary$pval), ]), n = 100))

# plot histogram of pvalues
plot(density(fit.afe$pval))  # NICE


# What about sitecounts? --------------------------------------------------

N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_matrix"
suffix <- "sitecounts.merged.matrix"
N <- LoadSitecountsEncodeAll(maindir = N.dir, suffix = suffix, with.ensemblid = FALSE)  # merged by gene


# Prepare sitecounts ------------------------------------------------------



# Let's look at Insig2 ----------------------------------------------------

jgene <- "Tars"
jgene <- "Galnt11"
jgene <- "Tcp1"
jgene <- "Ddc"
jgene <- "Trim2"
jgene <- "Insig2"
jmotif <- "RORA.p2"
# jmotif <- as.character(sample(x = unique(N$motif), size = 1))

N.sub <- subset(N, gene == jgene)
N.sub.motif <- subset(N, motif == "RORA.p2")
dat.rhyth.sub <- subset(fits.relamp, gene == jgene & tissue %in% tissues)
dat.rhyth.sub.motif <- subset(fits.relamp, tissue %in% tissues)
rownames(dat.rhyth.sub) <- dat.rhyth.sub$tissues  # indexing

# ggplot(N.sub, aes(x = tissue, y = motevo.value)) + geom_bar(stat = "identity") + ggtitle(paste(jgene, jmotif))

ampdic <- setNames(object = dat.rhyth.sub$relamp, nm = tissuesdic[as.character(dat.rhyth.sub$tissue)])
ampdic.motif <- setNames(object = dat.rhyth.sub.motif$relamp, nm = paste(tissuesdic[as.character(dat.rhyth.sub.motif$tissue)], dat.rhyth.sub.motif$gene, sep = ';'))
# add relamp info into N.sub
N.sub$relamp <- ampdic[as.character(N.sub$tissue)]
N.sub.motif$relamp <- ampdic.motif[paste(as.character(N.sub.motif$tissue), N.sub.motif$gene, sep = ';')]

# Quantify each motif by distance -----------------------------------------

fit.motif <- N.sub %>%
  group_by(motif) %>%
  do(FitMotifAmp(.))

fit.motif <- data.frame(fit.motif)
rownames(fit.motif) <- fit.motif$motif
head(data.frame(fit.motif[order(fit.motif$pval), ]), n = 50)
head(data.frame(fit.motif[order(abs(fit.motif$relamp), decreasing = TRUE), ]), n = 50)

ggplot(subset(N.sub, motif == "RORA.p2"), aes(y = motevo.value / sum(motevo.value), x = relamp)) + geom_point() + geom_smooth(method = "lm")


# Quantify a motif across all genes ---------------------------------------

N <- N %>%
  group_by(gene, tissue) %>%
  mutate(motevo.value.norm = motevo.value / sum(motevo.value))

rhythmic.genes <- as.character(head(data.frame(fit.afe.summary[order(fit.afe.summary$pval), ]), n = 100)$gene)
N.sub.motif <- subset(N, gene %in% rhythmic.genes & tissue %in% encodetissues)
dat.rhyth.sub.motif <- subset(fits.relamp, gene %in% rhythmic.genes & tissue %in% tissues)

N.sub.motif$relamp <- mapply(function(tissue, gene, dic) dic[[paste(tissue, gene, sep=";")]], 
                             as.character(N.sub.motif$tissue), as.character(N.sub.motif$gene),
                             MoreArgs = list(dic = ampdic.motif))

ggplot(subset(N.sub.motif, motif == "RORA.p2"), aes(y = motevo.value.norm, x = relamp)) + geom_point(alpha = 0.1) + geom_smooth(method = "lm")


# Quantify a motif for a specific gene ------------------------------------

genes <- c("Insig2", "Slc45a3", "Prkd3", "Gne")  # liver specific
genes <- c("Ddc", "Upp2", "Asl", "Prkd3", "Neu2")  # kidney specific
N.sub <- subset(N, gene %in% genes & tissue %in% encodetissues)
dat.rhyth.sub <- subset(fits.relamp, gene %in% genes & tissue %in% tissues)

N.sub$relamp <- mapply(function(tissue, gene, dic) dic[[paste(tissue, gene, sep=";")]],
                       as.character(N.sub$tissue), as.character(N.sub$gene),
                       MoreArgs = list(dic = ampdic.motif))

fit.motif <- N.sub %>%
  group_by(motif) %>%
  do(FitMotifAmp(.)) %>%
  arrange(pval)
(fit.motif[order(fit.motif$pval), ])
head(data.frame(fit.motif[order(fit.motif$relamp, decreasing = TRUE), ]), n = 100)
head(data.frame(fit.motif[order(fit.motif$pval, decreasing = FALSE), ]), n = 100)

ggplot(subset(N.sub, motif == "SREBF1.2.p2"), aes(x = motevo.value, y = relamp)) + geom_point() + geom_smooth(method = "lm")
