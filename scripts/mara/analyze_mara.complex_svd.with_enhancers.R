# Jake Yeung
# analyze_mara.complex_svd.with_enhancers.R (output from /home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/merged_dhs)
# 2015-09-10

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")

library(dplyr)  # problems with plyr old package
library(reshape2)

PenalizeNoise <- function(act.complex){
  # Create exprs.adj based on penalizing by noise
  act.complex$exprs.adj <- act.complex$exprs.transformed * act.complex$frac.weight
  act.complex$mod.exprs <- Mod(act.complex$exprs.transformed)
  act.complex$mod.exprs.adj <- Mod(act.complex$exprs.adj)
  return(act.complex)
}

# Load --------------------------------------------------------------------

indirs <- list(all="/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_dhs_multigene/expressed_genes_deseq_int.centeredTRUE",
               high="/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_dhs_multigene/entropy.high_entropy.centeredTRUE",
               low="/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_dhs_multigene/entropy.low_entropy.centeredTRUE")

act.longs <- lapply(indirs, LoadActivitiesLong)

# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "ATF2.p2")
lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "RORA.p2")
# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "HNF4A_NR2F1.2.p2")
lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "bHLH_family.p2")
# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "MYOD1.p2")
# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "ONECUT1.2.p2")
# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "HNF1A.p2")
# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "HNF4A_NR2F1.2.p2")

# To complex --------------------------------------------------------------

act.complexes <- lapply(act.longs, TemporalToFrequencyDatLong, period = 24, n = 8, interval = 6, add.entropy.method = "array")

act.complexes <- lapply(act.complexes, PenalizeNoise)

# act.complexes$noliver <- subset(act.complexes$noliver, ! tissue %in% c("Liver"))
# act.complexes$noliverbfatmus <- subset(act.complexes$noliverbfatmus, ! tissue %in% c("Liver", "BFAT", "Mus"))

# SVD ---------------------------------------------------------------------

s.acts <- lapply(act.complexes, SvdOnComplex, value.var = "exprs.adj")

lapply(s.acts, PlotFirstNComponents, comps = c(1, 2))


