# Jake Yeung
# analyze_mara.complex_svd.R
# 2015-09-09 
# Remove some tissues (Liver or Liver,BFAT,Mus) to see what happens to low and high entropy genes

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

# indirs <- list(all="/home/yeung/projects/tissue-specificity/results/MARA/entropy_nonadj_var_by_2_batch/entropy.low_entropy.centeredTRUE",
#                noliver="/home/yeung/projects/tissue-specificity/results/MARA/entropy_genes_split_by_2.no_liver_bfat_mus/entropy.low_entropy.centeredTRUE",
#                noliverbfatmus="/home/yeung/projects/tissue-specificity/results/MARA/entropy_genes_split_by_2.no_liver/entropy.low_entropy.centeredTRUE")
indirs <- list(all="/home/yeung/projects/tissue-specificity/results/MARA/entropy_nonadj_var_by_2_batch/entropy.high_entropy.centeredTRUE",
               noliver="/home/yeung/projects/tissue-specificity/results/MARA/entropy_genes_split_by_2.no_liver/entropy.high_entropy.centeredTRUE",
               noliverbfatmus="/home/yeung/projects/tissue-specificity/results/MARA/entropy_genes_split_by_2.no_liver_bfat_mus/entropy.high_entropy.centeredTRUE")
# indirs <- list(all="/home/yeung/projects/tissue-specificity/results/MARA/entropy_genes_split_by_3.all_tissues/entropy.high_entropy.centeredTRUE")
# indirs <- list(all="/home/yeung/projects/tissue-specificity/results/MARA/entropy_nonadj_var_by_2_batch/entropy.high_entropy.centeredTRUE")

act.longs <- lapply(indirs, LoadActivitiesLong)

# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "ATF2.p2")
lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "RORA.p2")
# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "HNF4A_NR2F1.2.p2")
# lapply(act.longs, function(dat, jgene) PlotActivitiesWithSE(subset(dat, gene == jgene)), jgene = "bHLH_family.p2")

# To complex --------------------------------------------------------------

act.complexes <- lapply(act.longs, TemporalToFrequencyDatLong, period = 24, n = 8, interval = 6, add.entropy.method = "array")

act.complexes <- lapply(act.complexes, PenalizeNoise)

# act.complexes$noliver <- subset(act.complexes$noliver, ! tissue %in% c("Liver"))
# act.complexes$noliverbfatmus <- subset(act.complexes$noliverbfatmus, ! tissue %in% c("Liver", "BFAT", "Mus"))

# SVD ---------------------------------------------------------------------

s.acts <- lapply(act.complexes, SvdOnComplex, value.var = "exprs.adj")

lapply(s.acts, PlotFirstNComponents, comps = c(1,2))


