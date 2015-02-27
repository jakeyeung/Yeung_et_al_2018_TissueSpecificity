# plot_exprs_from_list.R
# February 25 2015
library(ggplot2)
setwd("/home/yeung/projects/tissue-specificity")

scripts.dir <- "scripts"
funcs.dir <- "functions"
source(file.path(scripts.dir, funcs.dir, "LoadArrayRnaSeq.R"))
source(file.path(scripts.dir, funcs.dir, "PlotGeneAcrossTissues.R"))
source(file.path(scripts.dir, funcs.dir, "RemoveExtension.R"))
source(file.path(scripts.dir, funcs.dir, "ReadListToVector.R"))

dat <- LoadArrayRnaSeq()

fnames <- ReadListToVector("plots/nconds/7_conds_filtered_02_bicw/files.txt")

plot.dir <- "plots/nconds/7_conds_filtered_02_bicw/exprs_plots"

for (fname in fnames){
	jgenes <- ReadListToVector(fname)

	out.fname <- RemoveExtension(basename(fname))

	dat.sub <- subset(dat, gene %in% jgenes)
	
	pdf(file.path(plot.dir, paste0(out.fname, ".pdf")))
	for (jgene in jgenes){
	  print(jgene)
	  print(PlotGeneAcrossTissues(subset(dat.sub, gene == jgene), jtitle = jgene))
	}
	dev.off()
}
