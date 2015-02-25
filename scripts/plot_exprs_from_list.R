# plot_exprs_from_list.R
# February 25 2015

scripts.dir <- "scripts"
funcs.dir <- "functions"
source(file.path(scripts.dir, funcs.dir, "LoadArrayRnaSeq.R"))
source(file.path(scripts.dir, funcs.dir, "PlotGeneAcrossTissues.R"))

dat <- LoadArrayRnaSeq()

fname <- "/home/yeung/projects/tissue-specificity/plots/nconds/7_conds_filtered/7_conds_filtered2.txt"
jgenes <- ReadListToVector(fname)

dat.sub <- subset(dat, gene %in% jgenes)

for (jgene in jgenes){
  PlotGeneAcrossTissues(subset(dat.sub, gene == jgene), jtitle = jgene)
}
