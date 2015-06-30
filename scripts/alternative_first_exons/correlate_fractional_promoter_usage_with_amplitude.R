# Jake Yeung
# 2015-06-26
# Instead of cutoffs of "rhythmic" vs "not", we should just model fractional 
# promoter usage with relative amplitude

library(dplyr)
library(hash)
library(ggplot2)

source("scripts/functions/LoadKallisto.R")
source("scripts/functions/LoadLong.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/TissueSpecificRhythmicsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LoadArrayRnaSeq.R")

# Functions ---------------------------------------------------------------


# Load data ---------------------------------------------------------------

tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")
dat.long <- LoadArrayRnaSeq()

dat.rhyth <- FitRhythmicDatLong(dat.long)

dat.rhyth.relamp <- GetRelamp(fits = dat.rhyth, max.pval = 1e-3)

# hash for fast access
keys <- paste(dat.rhyth.relamp$tissue, dat.rhyth.relamp$gene, sep = ";")
relampdic <- hash(keys, dat.rhyth.relamp$relamp)
avgexprsdic <- hash(keys, dat.rhyth.relamp$int.rnaseq)


# Find tissue-specific genes ----------------------------------------------

pval.min <- 1e-5
pval.max <- 0.05
relamp.max <- 0.1
mean.cutoff <- 6

dat.rhyth.relamp$is.rhythmic <- apply(as.matrix(dat.rhyth.relamp), 1, IsRhythmicApply, 
                                      pval.min = pval.min, pval.max = pval.max, relamp.max = relamp.max, cutoff = mean.cutoff)

# Do not consider Cere, Hypo and BS into this -----------------------------

filter.tissues <- c()
filter.tissues <- c('Cere', 'BS', 'Hypo')
# 
dat.rhyth.relamp <- subset(dat.rhyth.relamp, !tissue %in% filter.tissues)


# Get tissue-specific genes -----------------------------------------------

# only tissue-specific if it contains TRUE and FALSE in a gene
dat.rhyth.relamp <- dat.rhyth.relamp %>%
  group_by(gene) %>%
  do(IsTissueSpecificDplyr(., include.betweens = FALSE))

tissue.spec <- unique(subset(dat.rhyth.relamp, is.tiss.spec == TRUE)$gene)
(length(tissue.spec))

# Get fractional usage ----------------------------------------------------

tpm.filt <- subset(tpm.merged, !tissue %in% filter.tissues)
# tpm.filt <- tpm.merged
tpm.afe <- tpm.filt %>%
  group_by(gene_name, tissue, time) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1))


# Model fractional isoform usage to rhythmic oscillations -----------------
tested.genes <- unique(dat.rhyth.relamp$gene)  # 20013 genes

tpm.avg <- tpm.afe %>%
  subset(., !tissue %in% filter.tissues & gene_name %in% tested.genes) %>%
  group_by(gene_name, transcript_id, tissue) %>%
  summarise(tpm_norm.avg = mean(tpm_normalized), tpm_norm.var = var(tpm_normalized))

tpm.avg$relamp <- mapply(function(tiss, gene) relampdic[[paste(tiss, gene, sep=";")]],
                         as.character(tpm.avg$tissue), 
                         as.character(tpm.avg$gene_name))

tpm.avg$int.rnaseq <- mapply(function(tiss, gene) avgexprsdic[[paste(tiss, gene, sep=";")]],
                             as.character(tpm.avg$tissue),
                             as.character(tpm.avg$gene_name))

tpm.avg.filt <- subset(tpm.avg, gene_name %in% tissue.spec & tpm_norm.var > 0 & int.rnaseq > 4)
tpm.fit <- tpm.avg.filt %>%
  group_by(gene_name, transcript_id) %>%
  do(FitPromoterUsageToAmplitude(.))

tpm.summary <- tpm.fit %>%
  group_by(gene_name) %>%
  do(SubsetMinPval(jdf = .))

tpm.summary$pval.adj <- p.adjust(tpm.summary$pval)

pval.cutoff <- 0.05
tpm.summary.filt <- subset(tpm.summary, pval <= pval.cutoff)
sig.hits <- tpm.summary.filt$gene_name

head(data.frame(tpm.summary[order(tpm.summary$pval), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$relamp.range, decreasing = TRUE), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$tpm_norm.range, decreasing = TRUE), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$tpm_norm.range, decreasing = TRUE), ]), n = 100)


tpm.summary.by.tpmrange <- data.frame(tpm.summary.filt[order(tpm.summary.filt$tpm_norm.range, decreasing = TRUE), which(colnames(tpm.summary.filt) != "transcript_id")])
(tpm.summary.by.tpmrange)

tpm.summary.nothits <- subset(tpm.summary, pval > pval.cutoff)
(head(data.frame(tpm.summary.nothits[order(tpm.summary.nothits$relamp.range, decreasing = TRUE), which(colnames(tpm.summary.nothits) != "transcript_id")]), n = 150))

# Calculate how many genes are alt prom and how many are sig hits ---------

tissue.spec.alt.prom <- unique(subset(tpm.avg, gene_name %in% tissue.spec & tpm_norm.var > 0 & int.rnaseq > 4)$gene_name)
sprintf("%s/%s genes correlate with alt promoter usage: %s",length(sig.hits), length(tissue.spec.alt.prom), length(sig.hits) / length(tissue.spec.alt.prom))

# Sanity checks -----------------------------------------------------------

# check hits
jgene <- "Insc"
jgene <- "Aqp3"
jgene <- "Slc43a1"
jgene <- "Ehd4"
jgene <- "Ddc"
jgene <- "Insig2"
jgene <- "Myh7"
jgene <- "Tnnt1"
jgene <- "Neb"
jgene <- "Tpm3"
jgene <- "Csrp3"
jgene <- "Atp5g1"
jgene <- "Insig2"
jgene <- "Zmym3"
jgene <- "Rapgef4"
jgene <- "Tmem216"
jgene <- "Ninl"
jgene <- "Cacna1h"
jgene <- "Ivd"
jgene <- "Bbs9"
jgene <- "Mpc1"
jgene <- "Stat5b"
jgene <- "Slc45a3"
jgene <- "Upp2"
jgene <- "Caprin1"
jgene <- "Atad3a"
jgene <- "Insc"
jgene <- "Ivd"
jgene <- "Bbs9"
jgene <- "Mpc1"

# check not hits
jgene <- "Csrp3"
jgene <- "Ncoa7"
jgene <- "Slc38a3"
jgene <- "Cdkn1a"
jgene <- "Usp13"

tpm.sub <- subset(tpm.fit, gene_name == jgene)
(tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ])
jtranscript <- tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1]
data.frame(subset(tpm.afe, transcript_id == jtranscript))
PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
ggplot(subset(tpm.avg.filt, transcript_id == jtranscript), aes(y = tpm_norm.avg, x = relamp, label = tissue)) + geom_point() + geom_text() + ggtitle(jgene) + geom_smooth(method = "lm")

