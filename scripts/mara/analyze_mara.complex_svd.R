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
act.long <- LoadActivitiesLong(indir)

PlotActivitiesWithSE(subset(act.long, gene == "RORA.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "LEF1_TCF7_TCF7L1.2.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "HNF1A.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "REST.p3"))
PlotActivitiesWithSE(subset(act.long, gene == "MEF2.A.B.C.D..p2" & tissue %in% c("Mus", "Liver")))
PlotActivitiesWithSE(subset(act.long, gene == "HNF4A_NR2F1.2.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "RXRG_dimer.p3"))
PlotActivitiesWithSE(subset(act.long, gene == "RUNX1..3.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "SRF.p3"))

# To complex --------------------------------------------------------------

act.complex <- TemporalToFrequencyDatLong(act.long, period = 24, n = 8, interval = 6, add.entropy.method = "array")

act.complex$exprs.adj <- act.complex$exprs.transformed * act.complex$frac.weight
act.complex$mod.exprs <- Mod(act.complex$exprs.transformed)
act.complex$mod.exprs.adj <- Mod(act.complex$exprs.adj)


# Plot examples -----------------------------------------------------------

jgene <- "KLF4.p3"
act.sub <- subset(act.complex, gene == "RUNX1..3.p2")
act.sub <- subset(act.complex, gene == jgene)
print(PlotComplex2(act.sub$exprs.transformed, labels = act.sub$tissue))
print(PlotComplex2(act.sub$exprs.adj, labels = act.sub$tissue))

# SVD ---------------------------------------------------------------------

s.act <- SvdOnComplex(act.complex, value.var = "exprs.adj")

for (comp in seq(3)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp)
  multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
}
