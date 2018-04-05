# 2016-06-25
# Jake Yeung

library(dplyr)
library(reshape2)
library(hash)
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/BiomartFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)


dat.livkid <- LoadLivKid()


# Remove bad genes --------------------------------------------------------

dat.long <- subset(dat.long, !is.na(gene))

# Remove pseudogenes
genes <- unique(as.character(dat.long$gene))
genes.biotype <- AnnotatePseudogenes(genes, return.original = FALSE)
genes.length <- AnnotateTranscriptLength(genes, return.original = FALSE)
biotype.hash <- hash(genes, genes.biotype)
length.hash <- hash(genes, genes.length)
dat.long$gbiotype <- sapply(as.character(dat.long$gene), function(g) biotype.hash[[g]])
dat.long$glength <- sapply(as.character(dat.long$gene), function(g) length.hash[[g]])
dat.long <- subset(dat.long, gbiotype == "protein_coding" & glength > 250)


dat.mean <- dat.long %>%
  group_by(gene) %>%
  summarise(exprs.max = quantile(exprs, probs = 0.9))

plot(density(dat.mean$exprs.max))
jcutoff <- 1
abline(v=jcutoff)

genes.cut <- as.character(subset(dat.mean, exprs.max <= jcutoff)$gene)

dat.long <- subset(dat.long, ! gene %in% genes.cut)

# PCA ---------------------------------------------------------------------

dat.sub <- StaggeredTimepointsLivKid(dat.long)

M <- dcast(subset(dat.sub, !is.na(gene)), formula = gene ~ tissue + geno + time, value.var = "exprs")
# M <- dcast(dat.livkid, formula = gene ~ tissue + time, value.var = "exprs")

rownames(M) <- M$gene
M$gene <- NULL

M <- t(scale(t(M), center = TRUE, scale = FALSE))

p <- prcomp(M, scale=FALSE, center=FALSE)

screeplot(p)

plot(p$rotation[, 1], p$rotation[, 2])
text(p$rotation[, 1], p$rotation[, 2], labels = rownames(p$rotation), cex = 1)

head(sort(abs(p$x[, 1]), decreasing=TRUE))

jgene <- "Gm14440"
jgene <- "Atf7ip"
jgene <- "Gm26602"
jgene <- "Ugt1a7c"
jgene <- "Tmem117"
PlotGeneTissuesWTKO(subset(dat.long, gene == jgene))
