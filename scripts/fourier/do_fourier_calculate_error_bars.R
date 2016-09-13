# 2016-09-09
# Jake Yeung
# do_fourier_calculate_error_bars.R

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/ActivitiesFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("/home/yeung/projects/posttranscriptional_regulation/functions/CosSineFunctions.R")

# Load --------------------------------------------------------------------

# get regulators: hogenesch 
outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.hogenesch"
outmain <- file.path(outbase, paste0("promoters.tissuewide.filteramp.0.15.mat"))
indir <- file.path(outmain, "expressed_genes_deseq_int.centeredTRUE")
act.long <- LoadActivitiesLong(indir, shorten.motif.name = TRUE)
# rename motifs based on the motifs with non-zero entries
omega <- 2 * pi / 24

PlotActivitiesWithSE(subset(act.long, gene == "bHLH_family"), jxlab = "CT") + facet_wrap(~tissue)

n <- 4
act.complex <- subset(act.long, tissue != "WFAT") %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE, propagate.errors = TRUE)) %>%
  mutate(amp = GetAmp(a = Re(exprs.transformed), b = Im(exprs.transformed), n = n),  # peak to trough
         phase = GetPhi(a = Re(exprs.transformed), b = Im(exprs.transformed), omega),  # in hours
         amp.se = GetAmp.se(a = Re(exprs.transformed), b = Im(exprs.transformed), sig.a = se.real, sig.b = se.im, n = n),
         phase.se = GetPhi.se(a = Re(exprs.transformed), b = Im(exprs.transformed), sig.a = se.real, sig.b = se.im, omega))

s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")
