# filter_genes.kallistoarray.R
# Jake Yeung
# 2015-06-24
# Filter kallisto long by a gene list, then write matrix to file


library(dplyr)
library(reshape2)

# Function ----------------------------------------------------------------

source("/home/yeung/projects/tissue-specificity/scripts/functions/WriteGeneListMat.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/LoadLong.R")

GetGenelist <- function(path){
  dat <- read.table(path, header = FALSE)
  return(as.character(dat$V1))
}


# Load --------------------------------------------------------------------

genelist.path <- "/home/yeung/projects/tissue-specificity/data/gene_lists/Kidney_rhythmic_genes.pvalmax0.05pvalmin1e-04relamp0.1.txt"
genelist.path <- "/home/yeung/projects/tissue-specificity/data/gene_lists/Lung_rhythmic_genes.pvalmax0.05pvalmin0.001relamp0.1.txt"

genelist <- GetGenelist(genelist.path)

# # ka.long
# array.path <- "/home/yeung/projects/tissue-specificity/data/exprs_matrices/array_adj_to_kallisto.slope07.txt"
# rna.seq.path <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
# # ka.long <- LoadLong(array.path, rna.seq.path, scale.factor = 100, pseudocount = 1)
# 
# load(file = "Robjs/dat.rhyth.kallisto.Robj")


PlotAmpPhase(dat = subset(dat.rhyth, pval <= 5e-2 & tissue == "Kidney" & gene %in% genelist))

