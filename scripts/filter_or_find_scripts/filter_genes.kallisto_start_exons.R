# 2015-07-03
# filter_genes.kallisto_start_exons.R

library(dplyr)

scale.factor <- 100
pseudocount <- 1

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadKallisto.R")
source("scripts/functions/LoadKallistoGene.R")
source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LoadArrayRnaSeq.R")


# Load data ---------------------------------------------------------------

tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")
# kallisto.wide <- LoadKallistoGene(inpath = "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt", form = "wide")
tpm.merged$experiment <- "rnaseq"

dat.long <- LoadArrayRnaSeq()

# Add up the genes --------------------------------------------------------

# I realized that you cannot find a "cutoff" if you look by start exons, therefore
# we will add up exprs of the gene in order to get a nice double expression
tpm.gene <- tpm.merged %>%
  group_by(gene_name, tissue, time) %>%
  summarize(tpm = sum(tpm))
head(tpm.gene)

tpm.gene$exprs <- log2(tpm.gene$tpm * scale.factor + pseudocount)

# find cutoff
exprs.log2 <- subset(tpm.gene, tpm > 0)$exprs
cutoff <- FindCutoff(x = exprs.log2, lambdas = c(0.2, 0.8), mus = c(1, 10), k = 2)  # max=5.572458


# Find expressed genes for each tissue ------------------------------------

tpm.mean <- tpm.gene %>%
  group_by(gene_name, tissue) %>%
  summarize(exprs.mean = mean(exprs)) %>%
  filter(exprs.mean > cutoff$maximum)


# Write to table ----------------------------------------------------------

WriteToFile <- function(dat, outdir){
  tissue <- as.character(dat$tissue)[1]
  outf <- file.path(outdir, paste0(tissue, ".txt"))
  sink(file = outf)
  for (gene in dat$gene_name){
    cat(gene)
    cat("\n")
  }
  sink()
  return(data.frame())
}

outdir="data/gene_lists/expressed_genes_by_tissue"
dir.create(outdir)
tpm.mean %>%
  group_by(tissue) %>%
  do(WriteToFile(., outdir=outdir))
