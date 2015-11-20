# 2015-11-10
# Jake Yeung
# cluster by promoters: sanity check by plotting

library(ggplot2)
library(dplyr)
library(hash)
library(ellipse)

setwd("/home/yeung/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/alternative_exon_usage/clustering"

# Functions ---------------------------------------------------------------

source("scripts/functions/PlotUCSC.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")


# load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)

load("Robjs/tpm.afe.avg.binary.Robj", verbose=T)
load("Robjs/tpm.gauss.bic_models.Robj", verbose=T)
# load("Robjs/tpm.fuzzy.bic_models.Robj", verbose=T)
# load("Robjs/tpm.gauss.filt.Robj", verbose=T)
# load("Robjs/tpm.gauss.filt.Robj", verbose=T)
# load("Robjs/tpm.gauss2.filt.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/tpm.merged.Robj", verbose=T)

# order -------------------------------------------------------------------

tpm.gauss <- tpm.gauss[order(tpm.gauss$center.dists, decreasing = T), ]
data.frame(head(subset(tpm.gauss, select = -sigs, n = 25)))
# Meis2, MyoD, Klf4


# Genome wide phenmonenon? ------------------------------------------------

plot(density(log(tpm.gauss$center.dists)))
plot(density(tpm.gauss$center.dists))

ggplot(tpm.gauss, aes(x = center.dists, y = intrascore2, label = gene_name)) + geom_point(alpha = 0.1) + scale_y_log10() + geom_text()

key <- as.character(fits.best$gene)
val <- fits.best$amp.avg
dic <- hash(key, val)
tpm.gauss$amp.avg <- sapply(tpm.gauss$gene_name, function(g) dic[[as.character(g)]])
ggplot(tpm.gauss, aes(x = center.dists, y = amp.avg, label = gene_name)) + geom_point(alpha = 0.1) + geom_text()

jgene <- "Npas2"
PromoterSpacePlots.nostics(subset(tpm.gauss, gene_name == jgene)$sigs[[1]], jgene = jgene)

tpm.hits <- subset(tpm.gauss, subset = center.dists > 0.2, select = -sigs)


# How many genes are rhythmic ---------------------------------------------

genes.rhyth <- subset(fits.best, model != "")$gene
genes.flat <- subset(fits.best, model == "")$gene


# how  many rhyth contain mulitple promoters -------------------------------------

# hash for rhyth or flat
key <- c(as.character(genes.rhyth), as.character(genes.flat))
val <- c(rep(TRUE, length(genes.rhyth)), rep(FALSE, length(genes.flat)))
rhyth.flat.dic <- hash(key, val)

# get number of promoters
tested.genes <- fits.best$gene
tpm.merged.sum <- subset(tpm.merged, tissue == "Adr" & time == 22 & gene_name %in% tested.genes) %>%
  group_by(gene_name) %>%
  summarise(n.proms = length(transcript_id))
tpm.merged.sum$is.rhyth <- sapply(tpm.merged.sum$gene_name, function(x) rhyth.flat.dic[[as.character(x)]])

ggplot(tpm.merged.sum,aes(x = n.proms)) + geom_histogram(binwidth=1)

tpm.merged.sum$has.multiprom <- sapply(tpm.merged.sum$n.proms, function(n){
  if (n > 1){
    return("MultiProm") 
  } else {
    return("Single")
  }
}
)

rhyth.multiprom.tbl <- table(tpm.merged.sum$is.rhyth, tpm.merged.sum$has.multiprom)

# fraction multiprom are hits
sum(rhyth.multiprom.tbl[, 1])
nrow(tpm.hits) / sum(rhyth.multiprom.tbl[, 1])  # ~1.3% not much

# what about genes with Aorta,BFAT?
genes.aortabfat <- subset(fits.best, model == "Aorta;BFAT")$gene
ggplot(subset(tpm.gauss, gene_name %in%genes.aortabfat), aes(x = center.dists, y = intrascore2, label = gene_name)) + geom_point(alpha = 0.1) + scale_y_log10() + geom_text()

genes.kidneyliver <- subset(fits.best, model == "Mus")$gene
ggplot(subset(tpm.gauss, gene_name %in%genes.kidneyliver), aes(x = center.dists, y = intrascore2, label = gene_name)) + geom_point(alpha = 0.1) + scale_y_log10() + geom_text()
