# plot_exprs_from_genelist.R
# March 3 2015
library(ggplot2)
setwd("/home/yeung/projects/tissue-specificity")

scripts.dir <- "scripts"
funcs.dir <- "functions"
source(file.path(scripts.dir, funcs.dir, "LoadArrayRnaSeq.R"))
source(file.path(scripts.dir, funcs.dir, "PlotGeneAcrossTissues.R"))
source(file.path(scripts.dir, funcs.dir, "RemoveExtension.R"))
source(file.path(scripts.dir, funcs.dir, "ReadListToVector.R"))

plot.dir <- "plots/component_1_gene_exprs.WFAT_removed.pdf"
genes <- ReadListToVector("results/component_1.WFAT.removed.txt")

dat <- LoadArrayRnaSeq()
dat.sub <- subset(dat, gene %in% genes)

pdf(plot.dir)

for (jgene in genes){
  print(jgene)
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == jgene), jtitle = jgene)) 
}
dev.off()
