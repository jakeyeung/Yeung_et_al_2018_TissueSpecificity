# analyze.gene.enrichment.R
# 
setwd("/home/yeung/projects/tissue-specificity")

# Functions ---------------------------------------------------------------
source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/FixGeneName.R")
source("scripts/functions/RemoveExtension.R")

# Loop through file names and get gene enrichment -------------------------

# fnames <- read.table("plots/nconds/7_conds_filtered_05_amp/files.txt")  # full file paths
# fnames <- as.character(unlist(fnames))  # loop-able
fnames <- ReadListToVector("plots/nconds/7_conds_filtered_05_amp/files.txt")
fnames <- unlist(sapply(fnames, FixGeneName))

# genes.bg <- read.tajble("plots/nconds/7_conds_filtered_05_amp/filtered_genes.txt")
# genes.bg <- as.character(unlist(genes.bg))
genes.bg <- ReadListToVector("plots/nconds/7_conds_filtered_05_amp/filtered_genes.txt")
genes.bg <- unlist(sapply(genes.bg, FixGeneName))

# precalculation of these objs makes things go faster
sym2entrez <- CreateSym2Entrez()
entrez2GO <- CreateEntrez2GO()

for (fname in fnames){
  genes.hit <- read.table(fname)
  genes.hit <- as.character(unlist(genes.hit))
  if (length(genes.hit) < 5){
    next
  }
  # Run ontologies for MF, BP and CC
  for (onto in c("BP", "MF", "CC")){
    fname.noext <- RemoveExtension(fname)
    fname.out <- paste0(fname.noext, '_', onto, ".GOtop")
    if (file.exists(fname.out)){
      next
    }
    res <- AnalyzeGeneEnrichment(genes.bg, genes.hit, 
                                 sym2entrez, entrez2GO, 
                                 convert.sym.to.entrez = TRUE, 
                                 which.ontology = onto, 
                                 write.path = fname.out, 
                                 node.size = 5, 
                                 FDR.cutoff = 0.05)
  }
}
