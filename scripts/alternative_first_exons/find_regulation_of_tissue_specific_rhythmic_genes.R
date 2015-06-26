# Find regulation of tissue specific rhythmic (TSR) genes 
# Bin TSR into categories: Liver, Kidney etc
# I do not regard the "alternative promoters" of the genes
# 2015-06-25

library(hash)
library(dplyr)
library(ggplot2)

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/LoadKallisto.R")
source("scripts/functions/LoadLong.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LoadSitecounts.R")
source("scripts/functions/FitMotifAmp.R")
source("scripts/functions/TissueSpecificRhythmicsFunctions.R")

# Functions ---------------------------------------------------------------



# Load data ---------------------------------------------------------------

tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")
dat.long <- LoadLong()

N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_matrix"
suffix <- "sitecounts.merged.matrix"
N <- LoadSitecountsEncodeAll(maindir = N.dir, suffix = suffix, with.ensemblid = FALSE, rename.tissues = TRUE)  # merged by gene

N <- N %>%
  group_by(gene, tissue) %>%
  mutate(motevo.value.norm = motevo.value / sum(motevo.value))

# Fit rhythmic ------------------------------------------------------------

dat.rhyth <- FitRhythmicDatLong(dat.long)

dat.rhyth.relamp <- GetRelamp(fits = dat.rhyth, max.pval = 1e-3)


# Find tissue-specific genes ----------------------------------------------

pval.min <- 1e-3
pval.max <- 0.05
relamp.max <- 0.1
mean.cutoff <- 6

dat.rhyth.relamp$is.rhythmic <- apply(as.matrix(dat.rhyth.relamp), 1, IsRhythmicApply, 
                                      pval.min = pval.min, pval.max = pval.max, relamp.max = relamp.max, cutoff = mean.cutoff)


# Do not consider Cere, Hypo and BS into this -----------------------------

filter.tissues <- c("Cere", "BS", "Hypo")

dat.rhyth.relamp <- subset(dat.rhyth.relamp, !tissue %in% filter.tissues)


# Get tissue-specific genes -----------------------------------------------

# only tissue-specific if it contains TRUE and FALSE in a gene
dat.rhyth.relamp <- dat.rhyth.relamp %>%
  group_by(gene) %>%
  do(IsTissueSpecificDplyr(.))

tissue.spec <- unique(subset(dat.rhyth.relamp, is.tiss.spec == TRUE)$gene)



# Write to file -----------------------------------------------------------

outfile <- paste0("~/projects/tissue-specificity/data/gene_lists/tissue_specific_rhythmic_genes.pvalmax", pval.max, "pvalmin", pval.min, "relamp", relamp.max, ".txt")
if (!file.exists(outfile)){
  sink(file = outfile)
  for (g in tissue.spec){
    cat(g)
    cat("\n")
  }
  sink()  
}

# Find Liver-only rhythmic genes ------------------------------------------

jtissue <- "Heart"
liver.rhyth.genes <- dat.rhyth.relamp %>%
  subset(., is.tiss.spec == TRUE) %>%
  group_by(gene) %>%
  do(GetTissueSpecific(., tissue = jtissue))

liver.rhyth.genes <- subset(liver.rhyth.genes, tiss.spec == FALSE)$gene
head(as.character(liver.rhyth.genes))

jgene <- sample(liver.rhyth.genes, 1)

PlotGeneAcrossTissues(subset(dat.long, gene == jgene))

# liver.rhyth.genes.df <- subset(dat.rhyth.relamp, is.tiss.spec == TRUE & tissue == jtissue & is.rhythmic == TRUE)
# liver.rhyth.genes.filt <- liver.rhyth.genes.df[order(liver.rhyth.genes.df$pval)[1:nrow(liver.rhyth.genes.df)], ]
# liver.rhyth.genes <- liver.rhyth.genes.filt$gene
# length(liver.rhyth.genes)

liver.outfile <- paste0("~/projects/tissue-specificity/data/gene_lists/", jtissue, "_rhythmic_genes.pvalmax", pval.max, "pvalmin", pval.min, "relamp", relamp.max, ".txt")
# if (!file.exists(liver.outfile)){
sink(file = liver.outfile)
for (g in liver.rhyth.genes){
  cat(g)
  cat("\n")
}
sink()  
# }


# Explore sitecounts ------------------------------------------------------


# Are rhythmic-only genes enriched for some factors? ----------------------

N.sub <- subset(N, gene %in% liver.rhyth.genes & !tissue %in% c("Cere"))
dat.rhyth.sub <- subset(dat.rhyth.relamp, gene %in% liver.rhyth.genes & !tissue %in% c("Cere"))

common.genes <- intersect(N.sub$gene, dat.rhyth.sub$gene)

keys <- paste(dat.rhyth.sub$tissue, dat.rhyth.sub$gene, sep = ";")
vals <- dat.rhyth.sub$relamp
relampdic <- hash(keys, vals)

N.sub$relamp <- mapply(function(tissue, gene, dic) dic[[paste(tissue, gene, sep=";")]],
                       N.sub$tissue, N.sub$gene,
                       MoreArgs = list(dic = relampdic))

# Fit relamp to motevo value

motif.fit <- N.sub %>%
  group_by(motif) %>%
  do(FitMotifAmp(.))

motif.fit[order(motif.fit$relamp, decreasing = TRUE), ]
motif.fit[order(motif.fit$pval), ]

# plot hits
jmotif <- "RORA.p2"

ggplot(subset(N.sub, motif == jmotif), aes(x = relamp, y = motevo.value.norm)) + geom_point(alpha = 0.1) + geom_smooth(method = "lm")
