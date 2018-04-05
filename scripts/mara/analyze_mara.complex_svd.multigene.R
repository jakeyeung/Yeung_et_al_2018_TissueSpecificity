# Jake Yeung
# analyze_mara.complex_svd.R
# 2016-04-03

rm(list=ls())

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
source("scripts/functions/PlotGeneAcrossTissues.R")
# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/RemoveP2Name.R")

library(dplyr)  # problems with plyr old package
library(reshape2)

# Load --------------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/results/MARA.multigene/expressed_genes_deseq_int.centeredTRUE.25000"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA.multigene.pairsonly/expressed_genes_deseq_int.centeredTRUE.25000"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long <- LoadActivitiesLong(indir)

PlotActivitiesWithSE(subset(act.long, gene == "CUX2.p2.RXRG_dimer.p3"))

# To complex --------------------------------------------------------------

act.complex <- TemporalToFrequencyDatLong(act.long, period = 24, n = 8, interval = 6, add.entropy.method = "array")

act.complex$exprs.adj <- act.complex$exprs.transformed * act.complex$frac.weight
act.complex$mod.exprs <- Mod(act.complex$exprs.transformed)
act.complex$mod.exprs.adj <- Mod(act.complex$exprs.adj)

# remove P2 
# act.complex$gene <- sapply(as.character(act.complex$gene), RemoveP2Name)

# act.complex <- subset(act.complex, !tissue %in% c("Liver"))

# Plot examples -----------------------------------------------------------

# jgene <- "HSF1.2.p2"
# jgene <- "RORA.p2"
# jgene <- "PAX5.p2"
# jgene <- "ONECUT1.2.p2"
# jgene <- "HIC1.p2"
# jgene <- "TLX1..3_NFIC.dimer..p2"
# jgene <- "FOSL2.p2"
# jgene <- "POU5F1.p2"
# # act.sub <- subset(act.complex, gene == "RUNX1..3.p2")
# act.sub <- subset(act.complex, gene == jgene)
# print(PlotComplex2(act.sub$exprs.transformed, labels = act.sub$tissue, title = jgene))

# SVD ---------------------------------------------------------------------

# s.act <- SvdOnComplex(subset(act.complex, ! gene %in% c("RORA.p2", "NFIL3.p2")), value.var = "exprs.adj")
s.act <- SvdOnComplex(act.complex, value.var = "exprs.adj")
# s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
for (comp in seq(1)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 14)
  multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
}

# list hits
hits.rhyth <- names(eigens.act$eigensamp[order(abs(Mod(eigens.act$eigensamp)), decreasing=TRUE)])


# Take top tissue-specific factors ----------------------------------------

act.mean <- subset(act.long, experiment == "rnaseq") %>%
  group_by(tissue, gene) %>%
  summarise(exprs.mean = mean(exprs))

act.diff <- act.mean %>%
  group_by(gene) %>%
  summarise(abs.diff = abs(diff(range(exprs.mean))))

act.diff <- act.diff[order(act.diff$abs.diff, decreasing = TRUE), ]

hits.tiss <- as.character(act.diff$gene)


# Take top hits -----------------------------------------------------------
# 
# top.rhyth <- hits.rhyth[1:20]
# top.tiss <- hits.tiss[1:20]

# sink(file = "/home/yeung/projects/tissue-specificity/data/gene_lists/motif_lists/top.rhyth.txt")
# for (m in top.rhyth){
#   cat(m); cat("\n")
# }
# sink()
# 
# sink(file = "/home/yeung/projects/tissue-specificity/data/gene_lists/motif_lists/top.tiss.txt")
# for (m in top.tiss){
#   cat(m); cat("\n")
# }
# sink()
# 
# sink(file = "/home/yeung/projects/tissue-specificity/data/gene_lists/motif_lists/all_motifs.txt")
# for (m in unique(act.long$gene)){
#   cat(m); cat("\n")
# }
# sink()
