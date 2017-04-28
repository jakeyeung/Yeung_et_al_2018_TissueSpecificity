# 2017-04-28
# Jake Yeung
# Can we identify rhythmic regulators by filtering liver DHSs then running MARA?

rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)

source("/home/yeung/projects/sleep_deprivation/scripts/functions/MaraDownstream.R")

source("scripts/functions/PlotGeneAcrossTissues.R")

# Load liver DHS ----------------------------------------------------------

jweight <- 0.8  # take all liver DHSs to all genes in Liver_SV129
jweight <- 0  # take all liver DHSs to all genes in Liver_SV129

flatampmax <- 0.1

promoters.only <- FALSE
all.genes <- FALSE
suffix <- paste0(".weight.", jweight, ".flatampmax.", flatampmax, ".promoters.", promoters.only, ".all_genes.", all.genes)

inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.200.weight.", jweight, ".MergePeaks.FALSE.nullmodel.JI.flatampmax.", flatampmax, ".withNmatallNmatfreqs.RemoveZeroCounts.Robj")
load(inf, v=T)

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)

liver.genes.all <- as.character(subset(fits.long.filt, model == "Liver_SV129")$gene)

if (all.genes){
  liver.genes <- liver.genes.all
} else {
  liver.peaks <- unique(as.character(subset(N.mat.all, model %in% "rhyth")$peak))
  liver.genes <- unique(as.character(subset(N.mat.all, model %in% "rhyth")$gene))
  print(paste("N peaks:, ", length(liver.peaks)))
  print(paste("N genes:, ", length(liver.genes)))
}


# Get gene expression over time and genotypes -----------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)


# Load MARA output --------------------------------------------------------

maraoutdir <- paste0("/home/yeung/data/tissue_specificity/mara_results/mara_outputs", suffix, "/center.TRUE", suffix, "/centered.TRUE/Liver_SV129_module_Liver_Spec_DHSs")
mara <- LoadMaraOutput(maraoutdir)

zscores <- mara$zscores
act.long <- mara$act.long
act.long$samp <- NULL

act.long$tissue <- sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[1]])
act.long$time <- as.numeric(sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[2]]))
act.long$geno <- as.character(sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[3]]))

PlotGeneTissuesWTKO(subset(act.long, gene == "RORA.p2"))
PlotGeneTissuesWTKO(subset(act.long, gene == "bHLH_family.p2"))
PlotGeneTissuesWTKO(subset(act.long, gene == "TFAP4.p2"))
PlotGeneTissuesWTKO(subset(act.long, gene == "HNF1A.p2"))
PlotGeneTissuesWTKO(subset(act.long, gene == "FOSL2.p2"))
PlotGeneTissuesWTKO(subset(act.long, gene == "AR.p2"))
PlotGeneTissuesWTKO(subset(act.long, gene == "FOXA2.p3"))


