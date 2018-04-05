# 2015-06-18
# find_oscillating_genes.kallistoarray.R
# We want to find oscillating genes, with p-value and amplitude proportional 
# to a core clock gene as a cutoff.

library(ggplot2)
library(dplyr)
# Functions ---------------------------------------------------------------

source("scripts/functions/PlotGeneTpm.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/LoadKallistoGene.R")
source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/WriteGeneListMat.R")

# Load matrix -------------------------------------------------------------

# array.path <- "data/exprs_matrices/array_adj_to_kallisto.slope07.txt"
# rna.seq.path <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
# ka.long <- LoadLong(array.path, rna.seq.path, scale.factor = 100, pseudocount = 1)
ka.long <- LoadArrayRnaSeq()

PlotGeneAcrossTissues(subset(ka.long, gene == "Arntl"))


# Center exprs for each gene across all samples ---------------------------

ka.centered <- ka.long %>%
  group_by(gene) %>%
  mutate(exprs.centered = scale(exprs, center = TRUE, scale = FALSE))

# Fit rhythmic ------------------------------------------------------------

fits <- FitRhythmicDatLong(ka.long)

save(fits, file = "Robjs/kallistoarray.fits.Robj")
# load(file = "Robjs/kallistoarray.fits.Robj")

# P-value adjuist ---------------------------------------------------------

fits$pval.adj <- NULL
fits <- fits %>%
  group_by(tissue) %>%
  mutate(pval.adj = p.adjust(pval)) %>%
  arrange(pval.adj)
# fits$pval.adj = p.adjust(fits$pval)

# Find most rhythmic gene: use this as relative amplitude for othe --------

max.pval <- 1e-3
# max.pvaladj <- 0.9

fits.mostrhythmic <- fits %>%
  group_by(tissue) %>%
  filter(pval <= max.pval) %>%
#   filter(pval.adj <= max.pvaladj) %>%
  do(FindMostRhythmic(.)) %>%
  data.frame(.)

rownames(fits.mostrhythmic) <- fits.mostrhythmic$tissue  # indexing


# Get amplitude relative to max amp ---------------------------------------

fits.relamp <- fits %>%
  group_by(tissue) %>%
  do(GetAmpRelToMax(., fits.mostrhythmic))

# List rhythmic genes -----------------------------------------------------

min.relamp <- 0.1  # fraction of largest amplitude of most rhythmic gene

fits.rhyth <- subset(fits.relamp, pval <= max.pval & amp >= min.relamp)
# fits.rhyth <- subset(fits.relamp, pval <= max.pval & relamp >= min.relamp)
# fits.rhyth <- subset(fits.relamp, pval.adj <= max.pvaladj & relamp >= min.relamp)

rhythmic.genes <- unique(fits.rhyth$gene)

print(length(rhythmic.genes))

# number of rhythmic genes in each tissue
rhyth.by.tiss <- fits.rhyth %>%
  group_by(tissue) %>%
  summarise(length(gene))

print(rhyth.by.tiss)


# Print gene list and expression matrix -----------------------------------

outdir <- paste0("/home/yeung/projects/tissue-specificity/data/gene_lists/rhythmic_genes_by_tissues_kallistoarray2", ".pval", max.pval, ".relamp", min.relamp)
dir.create(outdir)

fits.rhyth %>%
  group_by(tissue) %>%
  do(WriteGeneListMat(., ka.centered, outdir))
