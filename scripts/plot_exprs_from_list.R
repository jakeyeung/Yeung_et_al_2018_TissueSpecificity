# plot_exprs_from_list.R
# February 25 2015

scripts.dir <- "scripts"
funcs.dir <- "functions"
source(file.path(scripts.dir, funcs.dir, "LoadArrayRnaSeq.R"))
source(file.path(scripts.dir, funcs.dir, "PlotGeneAcrossTissues.R"))

dat <- LoadArrayRnaSeq()

fname <- "/home/yeung/projects/tissue-specificity/plots/nconds/7_conds_filtered/7_conds_filtered6.txt"
jgenes <- ReadListToVector(fname)

dat.sub <- subset(dat, gene %in% jgenes)

count <- 1
for (jgene in jgenes){
  print(jgene)
  print(PlotGeneAcrossTissues(subset(dat.sub, gene == jgene), jtitle = jgene))
  count <- count + 1
  if (count > 20){
    break
  }
}
