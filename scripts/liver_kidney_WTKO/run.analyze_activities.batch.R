# Jake Yeung
# run.analyze_activities.batch.R
# Analyze activities in large batches
# 2016-07-11

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

# dirmain <- "/home/yeung/projects/tissue-specificity/results/MARA/liver_kidney_wt_ko_cross_prods_analysis"
# dirmain <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_prods_tspeaks_TRUE"
# dirmain <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/cross_prods_tspeaks_cutofflow_05"
# dirmain <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.limitprod"
args <- commandArgs(trailingOnly=TRUE)
dirmain <- args[[1]]
outf <- args[[2]]

dirmainbase <- basename(dirmain)
fname <- paste0(dirmainbase, ".pdf")

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source('scripts/functions/LiverKidneyFunctions.R')

library(dplyr)  # problems with plyr old package
library(reshape2)
library(ggplot2)
library(ggrepel)

GetParam <- function(s, i){
  # expected string: "atger_kidney.dist_40000.cutoff_3.model_Liver_SV129.method_g=1001.cross_TRUE"
  # return Liver_SV129 if i = 4 for example
  param <- paste(strsplit(strsplit(s, "\\.")[[1]][[i]], "_")[[1]][-1], collapse = "_")
  return(param)
}

# Load --------------------------------------------------------------------


indirs <- list.dirs(dirmain, full.names = TRUE, recursive = FALSE)
print(indirs)

# outdir <- file.path("plots", "MARA", dirmainbase)

# print(paste("Outdir:", outdir))
# dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# pdf(file.path(outdir, fname))
pdf(outf)
for (indirmain in indirs){
  print(indirmain)
  indir <- list.dirs(indirmain, full.names = TRUE, recursive = FALSE)
  print(indir)
  act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = FALSE)

  omega <- 2 * pi / 24
  act.complex <- act.long %>%
    group_by(gene, tissue) %>%
    do(ProjectToFrequency2(., omega, add.tissue=TRUE))

  s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  jtitle <- gsub(pattern = "\\.", replacement = "\n", basename(indirmain))
  for (comp in seq(1)){
    eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 25, jtitle = jtitle)
    print(eigens.act$u.plot)
    # multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
  }
  # find RORA.p2.FOXA2.p3
  for (mtf in c("RORA", "SRF")){
    x <- Mod(eigens.act$eigensamp)
    # filter for RORA guys
    x <- x[grepl(mtf, names(x))]
    x <- x[order(x, decreasing = TRUE)]
    print(qplot(x = seq(length(x)), y = x, label = names(x), geom = "point") +  geom_text_repel() + theme_bw() + ggtitle(jtitle))
  }
}
dev.off()
