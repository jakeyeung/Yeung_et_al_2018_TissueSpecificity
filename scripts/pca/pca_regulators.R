# 2015-07-09
# Jake Yeung

library(reshape2)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/GetTFs.R")


# PCA on regulators -------------------------------------------------------

dat.long <- LoadArrayRnaSeq()

tfs <- GetTFs()
rbps <- ReadListToVector(fname = "/home/yeung/data/cisbp-rbp/RBP_names.txt", HEADER = TRUE)


# PCA on transcription factors --------------------------------------------

dat.tfs <- subset(dat.long, gene %in% tfs & experiment == "rnaseq")

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


# PCA on RBPs -------------------------------------------------------------

dat.rbps <- subset(dat.long, gene %in% rbps & experiment == "rnaseq")

dat.rbps.mean <- dat.rbps %>%
  group_by(gene, tissue) %>%
  summarise(Mean = mean(exprs))

M.rbps <- dcast(dat.rbps.mean, gene ~ tissue, value.var = "Mean")
rownames(M.rbps) <- M.rbps$gene
M.rbps$gene <- NULL

M.rbps <- t(scale(t(M.rbps), center = TRUE, scale = FALSE))
s.rbps <- prcomp(M.rbps, center = FALSE, scale. = FALSE)
screeplot(s.rbps)
biplot(s.rbps, choices = c(1, 2))

