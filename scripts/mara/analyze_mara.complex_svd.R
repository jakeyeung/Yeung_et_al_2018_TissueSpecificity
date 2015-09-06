# Jake Yeung
# analyze_mara.complex_svd.R
# 2015-07-17

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")

library(dplyr)  # problems with plyr old package
library(reshape2)

# Load --------------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_deseq_int"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/low_med_entropy_genes.centeredTRUE"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/low.entropy.genes.centeredTRUE/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/high.entropy.genes.centeredTRUE/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/lowmedhigh.entropy.genes.centeredTRUE/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/med.entropy.genes.centeredTRUE/"
act.long <- LoadActivitiesLong(indir)

PlotActivitiesWithSE(subset(act.long, gene == "RORA.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "NFIL3.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "bHLH_family.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "HIC1.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "TLX1..3_NFIC.dimer..p2"))
PlotActivitiesWithSE(subset(act.long, gene == "FOSL2.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "RFX1..5_RFXANK_RFXAP.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "MYOD1.p2"))
# # PlotActivitiesWithSE(subset(act.long, gene == "LEF1_TCF7_TCF7L1.2.p2"))
# # PlotActivitiesWithSE(subset(act.long, gene == "HNF1A.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "REST.p3"))
# PlotActivitiesWithSE(subset(act.long, gene == "MEF2.A.B.C.D..p2" & tissue %in% c("Mus", "BFAT")))
# PlotActivitiesWithSE(subset(act.long, gene == "HNF4A_NR2F1.2.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "RXRG_dimer.p3"))
# PlotActivitiesWithSE(subset(act.long, gene == "RUNX1..3.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "SRF.p3"))
# PlotActivitiesWithSE(subset(act.long, gene == "SREBF1.2.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "RFX1..5_RFXANK_RFXAP.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "NKX2.3_NKX2.5.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "ZNF423.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "EWSR1.FLI1.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "ETS1.2.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "ELK1.4_GABP.A.B1..p3"))
# PlotActivitiesWithSE(subset(act.long, gene == "PAX2.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "FOXA2.p3"))
# PlotActivitiesWithSE(subset(act.long, gene == "ONECUT1.2.p2"))

# To complex --------------------------------------------------------------

act.complex <- TemporalToFrequencyDatLong(act.long, period = 24, n = 8, interval = 6, add.entropy.method = "array")

act.complex$exprs.adj <- act.complex$exprs.transformed * act.complex$frac.weight
act.complex$mod.exprs <- Mod(act.complex$exprs.transformed)
act.complex$mod.exprs.adj <- Mod(act.complex$exprs.adj)

# act.complex <- subset(act.complex, !tissue %in% c("Liver"))

# Plot examples -----------------------------------------------------------

jgene <- "HSF1.2.p2"
jgene <- "RORA.p2"
jgene <- "PAX5.p2"
jgene <- "ONECUT1.2.p2"
jgene <- "HIC1.p2"
jgene <- "TLX1..3_NFIC.dimer..p2"
jgene <- "FOSL2.p2"
jgene <- "POU5F1.p2"
# act.sub <- subset(act.complex, gene == "RUNX1..3.p2")
act.sub <- subset(act.complex, gene == jgene)
print(PlotComplex2(act.sub$exprs.transformed, labels = act.sub$tissue, title = jgene))

# SVD ---------------------------------------------------------------------

s.act <- SvdOnComplex(act.complex, value.var = "exprs.adj")

for (comp in seq(3)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp)
  multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
}

