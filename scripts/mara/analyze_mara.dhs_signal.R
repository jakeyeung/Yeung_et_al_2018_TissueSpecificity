# Jake Yeung
# analyze_mara.complex_svd.R
# 2015-08-05

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotActivitiesFunctions.R")


# Load --------------------------------------------------------------------

# indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_DHS/MARA_dhs_gene_region_with_header"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_DHS/MARA_dhs_gene_region_cutoff_filter/"
# zscores.path <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_DHS/MARA_dhs_gene_region_with_header/Zscores"
zscores.path <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_DHS/MARA_dhs_gene_region_cutoff_filter/Zscores"
act.dhs.long <- LoadActivitiesLongDhs(indir, act.file = "Activities", se.file = "StandardError")
zscores <- read.table(zscores.path, header=FALSE, sep="\t", col.names = c("motif", "zscore"))
zscores <- zscores[order(zscores$zscore, decreasing = TRUE), ]

act.sum <- act.dhs.long %>%
            group_by(gene) %>%
            summarise(exprs.sum = sum(exprs))

# jgene <- "HNF4A_NR2F1.2.p2"
# jgene <- "REST.p3"
# jgene <- "RORA.p2"
# jgene <- "ONECUT1.2.p2"
# jgene <- "MEF2.A.B.C.D..p2"
# jgene <- "KLF4.p3"
# jgene <- "RORA.p2"
# jgene <- "RFX1..5_RFXANK_RFXAP.p2"
# jgene <- "FOXA2.p3"
# jgene <- "CTCF.p2"
# jgene <- "NR4A2.p2"
# jgene <- "GATA1..3.p2"
# jgene <- "TLX1..3_NFIC.dimer..p2"
# jgene <- "GATA6.p2"

pdf("plots/mara/mara_dhs_activities2.pdf")
for (m in zscores$motif){
  print(m)
  print(PlotActivitiesWithSE.dhs(subset(act.dhs.long, gene == m)))
}
dev.off()