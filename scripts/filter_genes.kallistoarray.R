# filter_genes.kallistoarray.R
# Jake Yeung
# 2015-06-24
# Filter kallisto long by a gene list, then write matrix to file

library(dplyr)
library(reshape2)

# Function ----------------------------------------------------------------

source("scripts/functions/WriteGeneListMat.R")

GetGenelist <- function(path){
  dat <- read.table(path, header = FALSE)
  return(as.character(dat$V1))
}


# Load --------------------------------------------------------------------

# gene list
genelist.path <- "/home/yeung/projects/tissue-specificity/data/gene_lists/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1/all.allgenelist"
genelist <- GetGenelist(genelist.path)

# ka.long
array.path <- "data/exprs_matrices/array_adj_to_kallisto.slope07.txt"
rna.seq.path <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
ka.long <- LoadLong(array.path, rna.seq.path, scale.factor = 100, pseudocount = 1)

# center expression
ka.centered <- ka.long %>%
  group_by(gene) %>%
  mutate(exprs.centered = scale(exprs, center = TRUE, scale = FALSE))


# write matrix expression
outdir <- "data/gene_lists/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1.union"

tissues <- unique(as.character(ka.centered$tissue))
for (tiss in tissues){
  fname.mat <- paste(tiss, "mat", sep=".")
  ka.sub <- subset(ka.centered, gene %in% genelist & tissue == tiss)
  m <- dcast(data = ka.sub, formula = gene ~ tissue + time + experiment, value.var = "exprs.centered")
  m <- FixColnamesForMARA(m)
  write.table(m, file = file.path(outdir, fname.mat), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
}