# 2015-07-09
# Jake Yeung

library(reshape2)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/GetTFs.R")


# PCA on regulators -------------------------------------------------------

dat.long <- LoadArrayRnaSeq(fix.rik.xgene = TRUE)

tfs.lst <- ReadListToVector(fname = "/home/yeung/projects/tissue-specificity/data/gene_lists/regulators/TFs_all.lst", HEADER = TRUE)
tfs <- unlist(sapply(tfs.lst, function(g) strsplit(g, split = ";")[[1]]), use.names = FALSE)

# PCA on transcription factors --------------------------------------------

filt.tiss <- c("Cere", "BS", "Hypo")
dat.tfs <- subset(dat.long, gene %in% tfs & experiment == "rnaseq" & !tissue %in% filt.tiss)

dat.tfs.mean <- dat.tfs %>%
  group_by(gene, tissue) %>%
  summarise(Mean = mean(exprs))

M.tfs <- dcast(dat.tfs.mean, gene ~ tissue, value.var = "Mean")
rownames(M.tfs) <- M.tfs$gene
M.tfs$gene <- NULL

M.tfs <- t(scale(t(M.tfs), center = TRUE, scale = FALSE))
s.tfs <- prcomp(M.tfs, center = FALSE, scale. = FALSE)
screeplot(s.tfs)
biplot(s.tfs, choices = c(1, 2))
biplot(s.tfs, choices = c(2, 3))
biplot(s.tfs, choices = c(3, 4))


# PCA on complex ----------------------------------------------------------

# load(file = "Robjs/dat.complex.maxexprs4.Robj")
s.tfs <- SvdOnComplex(subset(dat.complex, gene %in% tfs), value.var = "exprs.adj")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
for (comp in seq(3)){
  eigens.tfs <- GetEigens(s.tfs, period = 24, comp = comp)
  multiplot(eigens.tfs$v.plot, eigens.tfs$u.plot, layout = jlayout)
}

