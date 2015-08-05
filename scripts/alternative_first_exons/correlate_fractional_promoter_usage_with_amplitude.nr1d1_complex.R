# Jake Yeung
# 2015-06-26
# Normalize amplitude by Nr1d1 and also penalize noisy genes 

library(dplyr)
library(hash)
library(ggplot2)

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/LoadKallisto.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")


# Get data ----------------------------------------------------------------


dat.long <- LoadArrayRnaSeq()
dat.fit <- FitRhythmicDatLong(dat.long)
tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")

# Get complex matrix ------------------------------------------------------
ref.gene <- "Nr1d1"
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)
dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))

genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)
dat.complex <- TemporalToFrequencyDatLong(subset(dat.long, gene %in% genes.exprs), period = 24, n = 8, interval = 6, add.entropy.method = "array")

# adjust for "noise"
dat.complex$exprs.adj <- dat.complex$exprs.transformed * dat.complex$frac.weight

# normalize to reference gene (Nr1d1)
ref.amps <- subset(dat.complex, gene == ref.gene, select = c(tissue, gene, exprs.adj))
ref.amps.dic <- hash(ref.amps$tissue, ref.amps$exprs.adj)

dat.complex$exprs.adj.norm <- mapply(function(tiss, exprs) exprs / ref.amps.dic[[tiss]], as.character(dat.complex$tissue), dat.complex$exprs.adj)
dat.complex$mod.exprs.adj <- Mod(dat.complex$exprs.adj)
dat.complex$mod.exprs.adj.norm <- Mod(dat.complex$exprs.adj.norm)


# Annotate complex mat with mean value ------------------------------------

means.dic.keys <- paste(dat.fit.relamp$tissue, dat.fit.relamp$gene, sep = ";")
means.dic.vals <- dat.fit.relamp$int.rnaseq
means.dic <- hash(means.dic.keys, means.dic.vals)

dat.complex$mean <- mapply(function(gene, tissue) means.dic[[paste(c(tissue, gene), collapse = ";")]], 
                           as.character(dat.complex$gene), as.character(dat.complex$tissue))

# Get only genes where at least one tissue is rhythmic --------------------

rhyth.genes <- unique(subset(dat.complex, mod.exprs.adj.norm >= 0.1)$gene)


# Calculate frac usage ----------------------------------------------------


tpm.afe <- tpm.merged %>%
  filter(gene_name %in% rhyth.genes) %>%
  group_by(gene_name, tissue, time) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1))

tpm.avg <- tpm.afe %>%
  group_by(gene_name, transcript_id, tissue) %>%
  summarise(tpm_norm.avg = mean(tpm_normalized), tpm_norm.var = var(tpm_normalized)) %>%
  filter(tpm_norm.avg < 1 & tpm_norm.var > 0)


# Annotate rel amp and mean exprs -----------------------------------------


# annotate with normalized relative amplitude and mean expression
relamp.dic.keys <- paste(dat.complex$tissue, dat.complex$gene, sep = ";")
relamp.dic.vals <- dat.complex$mod.exprs.adj.norm
relamp.dic <- hash(relamp.dic.keys, relamp.dic.vals)

tpm.avg$relamp <- mapply(function(gene, tissue) relamp.dic[[paste(c(tissue, gene), collapse = ";")]], 
                              as.character(tpm.avg$gene_name), as.character(tpm.avg$tissue))
tpm.avg$int.rnaseq <- mapply(function(gene, tissue) means.dic[[paste(c(tissue, gene), collapse = ";")]],
                       as.character(tpm.avg$gene_name), as.character(tpm.avg$tissue))

# Model frac usage to amp -------------------------------------------------

tpm.avg.filt <- subset(tpm.avg, int.rnaseq > 4)
tpm.fit <- tpm.avg.filt %>%
  group_by(gene_name, transcript_id) %>%
  do(FitPromoterUsageToAmplitude(.))


# Show hits ---------------------------------------------------------------

tpm.summary <- tpm.fit %>%
  group_by(gene_name) %>%
  do(SubsetMinPval(jdf = .))

pval.cutoff <- 0.05
tpm.summary.filt <- subset(tpm.summary, pval <= pval.cutoff)
sig.hits <- tpm.summary.filt$gene_name

head(data.frame(tpm.summary[order(tpm.summary$pval), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$relamp.range, decreasing = TRUE), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$tpm_norm.range, decreasing = TRUE), which(colnames(tpm.summary) != "transcript_id")]), n = 100)
head(data.frame(tpm.summary[order(tpm.summary$tpm_norm.range, decreasing = TRUE), ]), n = 100)


# Hits --------------------------------------------------------------------

plot(density(tpm.summary$pval[which(!is.na(tpm.summary$pval))]))

hits <- head(data.frame(tpm.summary[order(tpm.summary$pval), which(colnames(tpm.summary) != "transcript_id")]), n = 200)$gene_name

pdf("plots/alternative_exon_usage/linear_model_complex_mat_top_hits.pdf")
tpm.hits <- subset(tpm.afe, gene_name %in% hits)
tpm.avg.filt.hits <- subset(tpm.avg.filt, gene_name %in% hits)
dat.long.hits <- subset(dat.long, gene %in% hits)
for (jgene in hits){
  tpm.sub <- subset(tpm.fit.hits, gene_name == jgene)
  jtranscript <- tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1]
  PlotDiagnostics2(subset(tpm.hits, gene_name == jgene), 
                   subset(tpm.avg.filt.hits, transcript_id == jtranscript), 
                   subset(dat.long.hits, gene == jgene), jgene, jtranscript)
}
dev.off()

jgene <- "Upp2"
jgene <- "Ddc"
jgene <- "Insig2"
jgene <- "Tpm3"
jgene <- "Prkd3"
jgene <- "Cstad"
jgene <- "Npas2"

tpm.sub <- subset(tpm.fit, gene_name == jgene)
# (tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ])
jtranscript <- tpm.sub[order(abs(tpm.sub$relamp), decreasing = TRUE), ]$transcript_id[1]

PlotDiagnostics2(tpm.afe, tpm.fit, dat.long, jgene, jtranscript)

# data.frame(subset(tpm.afe, transcript_id == jtranscript))
subset(tpm.summary, gene_name == jgene)
PlotGeneAcrossTissues(subset(dat.long, gene == jgene))
ggplot(subset(tpm.avg.filt, transcript_id == jtranscript), aes(y = tpm_norm.avg, x = relamp, label = tissue)) + geom_point() + geom_text() + ggtitle(jgene) + geom_smooth(method = "lm")
