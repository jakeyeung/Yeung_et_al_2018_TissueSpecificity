# 2015-06-18
# find_oscillating_genes.kallistoarray.R
# We want to find oscillating genes, with p-value and amplitude proportional 
# to a core clock gene as a cutoff.

library(ggplot2)

# Functions ---------------------------------------------------------------

source("scripts/functions/PlotGeneTpm.R")
source("scripts/functions/LoadArray.R")
source("scripts/functions/LoadKallistoGene.R")
source("scripts/functions/GetTissueTimes.R")

LoadLong <- function(array.path, rna.seq.path, pseudocount = 1e-5){
  kallisto.wide <- LoadKallistoGene(rna.seq.path, form = "wide")  # adjusts colnames to match array
  array.wide <- read.table(array.path)
  common.genes <- intersect(rownames(kallisto.wide), rownames(array.wide))
  kallisto.sub <- kallisto.wide[common.genes, ]
  array.sub <- array.wide[common.genes, ]
  
  tissues.rnaseq <- GetTissues(colnames(kallisto.sub), get_unique = FALSE)
  times.rnaseq <- GetTimes(colnames(kallisto.sub), get_unique = FALSE)
  tissues.array <- GetTissues(colnames(array.sub), get_unique = FALSE)
  times.array <- GetTimes(colnames(array.sub), get_unique = FALSE)
  
  ka.long <- data.frame(gene = c(rep(rownames(kallisto.sub), ncol(kallisto.sub)), rep(rownames(array.sub), ncol(array.sub))),
                        tissue = c(rep(tissues.rnaseq, each = nrow(kallisto.sub)), rep(tissues.array, each = nrow(array.sub))),
                        time = as.numeric(c(rep(times.rnaseq, each = nrow(kallisto.sub)), rep(times.array, each = nrow(array.sub)))),
                        tpm = c(unlist(kallisto.sub), unlist(array.sub)),
                        experiment = c(rep("rnaseq", nrow(kallisto.sub) * ncol(kallisto.sub)), rep("array", nrow(array.sub) * ncol(array.sub))))
  ka.long$tpm <- log2(ka.long$tpm + pseudocount)
  return(ka.long)
}

PlotGeneAcrossTissues2 <- function(dat, jgene){
  dat.sub <- subset(dat, gene == jgene)
  m <- ggplot(dat.sub, aes(x = time, y = tpm, colour = experiment, group = experiment, fill = experiment)) +
    geom_point() + geom_line() + facet_wrap(~tissue)
    ggtitle(jgene)
  print(m)
}

# Load matrix -------------------------------------------------------------

array.path <- "data/exprs_matrices/array_adj_to_kallisto.slope07.txt"
rna.seq.path <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
ka.long <- LoadLong(array.path, rna.seq.path)

PlotGeneAcrossTissues2(ka.long, "Arntl")
