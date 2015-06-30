# filter_genes.kallistoarray.R
# Jake Yeung
# 2015-06-24
# Filter kallisto long by a gene list, then write matrix to file


library(dplyr)
library(reshape2)

# Function ----------------------------------------------------------------

source("/home/yeung/projects/tissue-specificity/scripts/functions/WriteGeneListMat.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/LoadLong.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/LoadArrayRnaSeq.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/GrepRikGenes.R")

GetGenelist <- function(path){
  dat <- read.table(path, header = FALSE)
  return(as.character(dat$V1))
}


# Load --------------------------------------------------------------------

jargs <- commandArgs(trailingOnly = TRUE)
genelist.path <- jargs[1]
# genelist.path <- "/home/yeung/projects/tissue-specificity/data/gene_lists/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1/all.allgenelist"
outdir <- jargs[2]
normalize <- jargs[3]
if (is.na(normalize)){
  normalize <- FALSE
} else if (normalize == "TRUE"){
  normalize <- TRUE
} else if (normalize == "FALSE"){
  normalize <- FALSE
} else {
  normalize <- TRUE
  print('Defaulting normalize to TRUE')
}
print(paste("Centering data:", normalize))
# outdir <- "data/gene_lists/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1.union"
dir.create(outdir)

genelist <- GetGenelist(genelist.path)

# ka.long
ka.long <- LoadArrayRnaSeq()
ka.long$gene <- as.character(ka.long$gene)
Xgenes <- GrepRikGenes(ka.long$gene)
Xgenes.i <- which(ka.long$gene %in% Xgenes)
Xgenes.fixed <- RemoveX(as.character(Xgenes))
ka.long$gene[Xgenes.i] <- Xgenes.fixed
ka.long$gene <- as.factor(ka.long$gene)

if (normalize){
  # center expression
  ka.centered <- ka.long %>%
    group_by(gene) %>%
    mutate(exprs.centered = scale(exprs, center = TRUE, scale = FALSE))
} else {
  ka.centered <- ka.long
  ka.centered$exprs.centered = ka.long$exprs
  print("Did not center")
}


# write matrix expression

tissues <- unique(as.character(ka.centered$tissue))
for (tiss in tissues){
  fname.mat <- paste(tiss, "mat", sep=".")
  ka.sub <- subset(ka.centered, gene %in% genelist & tissue == tiss)
  m <- dcast(data = ka.sub, formula = gene ~ tissue + time + experiment, value.var = "exprs.centered")
  m <- FixColnamesForMARA(m)
  write.table(m, file = file.path(outdir, fname.mat), quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
}