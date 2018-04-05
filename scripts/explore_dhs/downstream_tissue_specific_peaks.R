# 2016-05-09
# Downstream of calculate_tissue(liver)_specific_peaks_by_genelist.R

rm(list=ls())

library(ggplot2)

source("scripts/functions/PlotFunctions.R")

# Load data ---------------------------------------------------------------

maindir="/home/yeung/projects/tissue-specificity/Robjs/dhs_peaks"
df.out.lst.merged <- lapply(list.files(maindir), function(fname){
  fpath <- file.path(maindir, fname)
  tiss <- strsplit(fname, "\\.")[[1]][[10]]
  print(fpath)
  print(tiss)
  load(fpath)
  df.out.lst.merged$tissue <- tiss
  return(df.out.lst.merged)
})
df.out.lst.merged <- do.call(rbind, df.out.lst.merged)
head(df.out.lst.merged)
jsub <- subset(df.out.lst.merged, gene.type %in% c("jgenes", "jgenes.flat.filt"))
newlab <- hash(as.character(subset(jsub, gene.type == "jgenes")$tissue), paste0(subset(jsub, gene.type == "jgenes")$tissue, "\nN=", subset(jsub, gene.type == "jgenes")$total.genes))
jsub$Tissue <- sapply(as.character(jsub$tissue), function(tiss) newlab[[tiss]])
jsub <- OrderDecreasing(jsub, "Tissue", "frac.n.spec.by.gene")
ggplot(jsub, aes(x = Tissue, y = frac.n.spec.by.gene, fill = gene.type)) + geom_bar(stat = "identity", position = "dodge") + theme_bw(24) + xlab("") + ylab("Fraction of genes with tissue-specific DHS within 5kb")

