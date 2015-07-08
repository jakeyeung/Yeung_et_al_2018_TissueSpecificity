# 2015-07-06
# Some bug with DHS peaks (need to filter low values)
# for now, let's focus on Tars 

setwd("~/projects/tissue-specificity/")

library(dplyr)
library(ggplot2)
library(reshape2)
library(PMA)

source("scripts/functions/LoadSitecounts.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/DifferentialSitecountsFunctions.R")

# Load --------------------------------------------------------------------

# N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo_tars"
N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo_tars_nokidneypeak"
suffix <- "filtered.100000.noN.tars.nokid.mat"
N <- LoadSitecountsEncodeAll(maindir = N.dir, tissues = c("Liver", "Kidney", "Cere", "Lung", "Heart", "Mus"),
                             suffix = suffix, with.ensemblid = FALSE, rename.tissues = FALSE)  # merged by gene
dat.long <- LoadArrayRnaSeq()
load(file = "Robjs/dat.rhyth.relamp.pvalmin1e-5.pvalmax0.05.relampmax.0.1.meancutoff6.Robj", verbose = T)

N <- N %>%
  group_by(gene, tissue) %>%
  mutate(motevo.value.norm = motevo.value / sum(motevo.value))

dhs.tiss <- unique(N$tissue)

tfs <- GetTFs()


# Init exprs gene ---------------------------------------------------------

X.exprs <- dcast(data = subset(dat.rhyth.relamp, gene %in% tfs & tissue %in% dhs.tiss), formula = gene ~ tissue, value.var = "int.rnaseq")
rownames(X.exprs) <- X.exprs$gene
X.exprs$gene <- NULL

# Only tars ---------------------------------------------------------------

jgene <- "Tars"

subset(N, gene == jgene & motif == "RORA.p2")
subset(N, gene == jgene & motif == "HNF1A.p2")
subset(N, gene == jgene & motif == "HIC1.p2")
subset(N, gene == jgene & motif == "AHR_ARNT_ARNT2.p2")

N.sub <- subset(N, gene == jgene)
X.motif <- dcast(data = N.sub, formula = motif ~ tissue, value.var = "motevo.value")
rownames(X.motif) <- X.motif$motif
X.motif$motif <- NULL
# replace NA with 0
X.motif[is.na(X.motif)] <- 0

# center stuff
jscale <- FALSE
jcenter <- TRUE
X.exprs.scaled <- ScaleRemoveInfs(X.exprs)
X.motif.scaled <- ScaleRemoveInfs(X.motif)

p.motif <- prcomp(X.motif.scaled, center = TRUE, scale. = FALSE)
biplot(p.motif, main = paste(jgene, "Motif PCA"), cex = c(0.5, 1.6), pch = 20) 


# Rotate by LIVER vector --------------------------------------------------

V <- p.motif$rotation[, c(1, 2)]
U <- p.motif$x[, c(1, 2)]

liver.vec <- V["Liver", ]
liver.vec.norm <- liver.vec / sqrt(sum(liver.vec ^ 2))

V.livproj <- V %*% liver.vec.norm
U.livproj <- U %*% liver.vec.norm
U.livproj <- U.livproj[order(U.livproj, decreasing = TRUE), ]
U.livproj <- U.livproj[which(abs(U.livproj) > 1)]

par(mar=c(10.1, 4.1, 4.1, 2.1))
barplot(U.livproj, names.arg = names(U.livproj), las = 2)


# Do CCA with penalties ---------------------------------------------------

jscale <- FALSE
jcenter <- TRUE
X.exprs.scaled <- ScaleRemoveInfs(X.exprs)
X.motif.scaled <- ScaleRemoveInfs(X.motif)

X.motif.exprs <- MatchColumns(X.motif.scaled, X.exprs.scaled)
perm.out <- CCA.permute(t(X.motif.exprs$X.motif), t(X.motif.exprs$X.exprs), typex="standard", typez="standard", standardize = F)
penaltyx <- perm.out$bestpenaltyx
penaltyz <- perm.out$bestpenaltyz
# penaltyx <- 1
# penaltyz <- 1
cca.out <- CCA(t(X.motif.exprs$X.motif), t(X.motif.exprs$X.exprs), typex="standard", typez="standard", K=2, penaltyx=penaltyx, penaltyz=penaltyz, standardize = F)
rownames(cca.out$u) <- rownames(X.motif.exprs$X.motif)
rownames(cca.out$v) <- rownames(X.motif.exprs$X.exprs)

# visualize with biplots
Xu.motif <- t(X.motif.exprs$X.motif) %*% cca.out$u
Xv.exprs <- t(X.motif.exprs$X.exprs) %*% cca.out$v
cor(Xu.motif, Xv.exprs)
biplot(cca.out$u[, 1:2], Xu.motif[, 1:2])
biplot(cca.out$v[, 1:2], Xv.exprs[, 1:2])


# Visualize component 1 by bar --------------------------------------------

par(mar=c(10.1, 4.1, 4.1, 2.1))
u.sorted <- cca.out$u[order(cca.out$u[, 1]), 1]
u.sorted <- u.sorted[which(abs(u.sorted) > 0.05)]
barplot(u.sorted, names.arg = names(u.sorted), las = 2)


