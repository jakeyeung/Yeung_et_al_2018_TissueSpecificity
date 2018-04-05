# 2015-11-10
# Jake Yeung
# cluster by promoters: sanity check by plotting

library(dplyr)
library(hash)
library(ellipse)

setwd("/home/yeung/projects/tissue-specificity")
outdir <- "/home/yeung/projects/tissue-specificity/plots/alternative_exon_usage/clustering"

# Functions ---------------------------------------------------------------

source("scripts/functions/PlotUCSC.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")

GetBedAcrossPromoters <- function(dat){
  chromos <- unique(dat$chromo)
  if (length(chromos == 1)){
    chromo = chromos[[1]]
  } else {
    warning("Duplicate chromosomes found..")
  }
  start.min <- min(dat$start)
  end.max <- max(dat$end)
  return(data.frame(chromo = chromo, start = start.min, end = end.max))
}

# load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)

load("Robjs/tpm.afe.avg.binary.Robj", verbose=T)
load("Robjs/tpm.gauss.bic_models.Robj", verbose=T)
# load("Robjs/tpm.fuzzy.bic_models.Robj", verbose=T)
# load("Robjs/tpm.gauss.filt.Robj", verbose=T)
# load("Robjs/tpm.gauss.filt.Robj", verbose=T)
# load("Robjs/tpm.gauss2.filt.Robj", verbose=T)
load("Robjs/tpm.merged.Robj", verbose=T)
load("Robjs/tpm.max.dists.all.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

# order -------------------------------------------------------------------

tpm.gauss <- tpm.gauss[order(tpm.gauss$center.dists, decreasing = T), ]
data.frame(head(subset(tpm.gauss, select = -sigs, n = 25)))
# Meis2, MyoD, Klf4


# plot hits: gauss BIC models ----------------------------------------------

top <- 500
tpm.hits <- tpm.gauss[1:top, ]
hits <- tpm.gauss$gene_name[1:top]

# sort by top distances
outname <- paste0("diagnostics_alt_prom_usage_genomewide_", top, ".pdf")
hits.genomewide <- as.character(max.dists.all$gene_name)
pdf(file.path(outdir, outname))
for (jgene in hits.genomewide){  
  print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene)))
  jsub <- subset(tpm.gauss, gene_name == jgene)
  if (nrow(jsub) > 0){
    PromoterSpacePlots.nostics(subset(tpm.gauss, gene_name == jgene)$sigs[[1]], jgene, draw.ellipse = T)
  } else {
    PlotPromoter(subset(tpm.afe.avg, gene_name == jgene), jtitle = jgene)
  }
}
dev.off()

outname <- paste0("diagnostics_", top, ".pdf")
pdf(file.path(outdir, outname))
for (jgene in hits){
  print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene))); PromoterSpacePlots.nostics(subset(tpm.gauss, gene_name == jgene)$sigs[[1]], jgene, draw.ellipse = T)
}
dev.off()

tpm.sub <- subset(tpm.merged, gene_name %in% hits, select = c(chromo, start, end, gene_name, tissue, time)) %>%
  group_by(gene_name) %>%
  do(GetBedAcrossPromoters(.))

keys <- as.character(tpm.hits$gene_name)
vals <- as.numeric(tpm.hits$center.dists)
dic <- hash(keys, vals)
# reorder by hits
tpm.sub$center.dists <- sapply(as.character(tpm.sub$gene_name), function(jgene) dic[[jgene]])
tpm.sub <- tpm.sub[order(tpm.sub$center.dists, decreasing = T), ]

print(tpm.sub)
outname2 <- paste0("hitsucsc_", top, ".pdf")
bedToUCSC(toPlot = tpm.sub, outpdf = file.path(outdir, outname2))


