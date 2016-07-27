# 2016-07-27
# Jake Yeung
# Downstream analysis liver and kidney 

rm(list=ls())

library(ggplot2)
library(dplyr)

source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
dat.orig <- dat.long
dat.long <- StaggeredTimepointsLivKid(dat.long)

jmod <- "all"
# outmain <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.Kidney_SV129,Kidney_BmalKO.g=1001"
outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod, ".g=1001")
indir <- file.path(outmain, "atger_with_kidney.bugfixed")
source("scripts/functions/LoadActivitiesLong.R")
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = FALSE)


# Get mean ----------------------------------------------------------------


act.long$tiss <- as.factor(sapply(as.character(act.long$tissue), function(jtiss) strsplit(jtiss, "_")[[1]][[1]]))
# rename motifs based on the motifs with non-zero entries

act.mean <- act.long %>%
  group_by(gene, geno, tiss) %>%
  summarise(exprs = mean(exprs))


# find largest tissue-specific motifs
act.tissmean <- act.mean %>%
  group_by(gene, tiss) %>%
  summarise(exprs = mean(exprs))

act.tissdiff <- act.tissmean %>%
  group_by(gene) %>%
  summarise(exprs.diff = diff(exprs)) %>%
  # arrange(desc(abs(exprs.diff)))
  arrange(exprs.diff)

print(head(act.tissdiff))  # MYBL2.p2 etc
jmotifs <- head(as.character(act.tissdiff$gene), n = 20)
plotdir <- "plots/mara_liver_kidney_modules"
pdf(file.path(plotdir, paste0("tissue_spec_motifs.", jmod, ".pdf")))
# jmotif <- "HSF1.2.p2"
for (jmotif in jmotifs){
  print(PlotActivitiesWithSE.rnaseq(subset(act.long, gene == jmotif)))
  tfs <- GetTFs(get.mat.only = TRUE)
  for (g in GetGenesFromMotifs(jmotif, tfs)){
    print(g)
    jsub <- subset(dat.orig, gene == g)
    if (nrow(jsub)) 
    print(PlotGeneTissuesWTKO(jsub, jtitle = g))
  }
}
dev.off()


# Find genes highly expressed in Kidney but not in Liver ------------------

dat.diff <- dat.long %>%
  group_by(tissue, gene) %>%
  summarise(exprs = mean(exprs)) %>%
  group_by(gene) %>%
  summarise(exprs.diff = diff(exprs)) %>%
  arrange(desc(exprs.diff))


source("scripts/functions/ReadListToVector.R")
tfs.lst <- GetTFs()
rbps <- ReadListToVector(fname = "/home/yeung/data/cisbp-rbp/RBP_names.txt", HEADER = TRUE)

subset(dat.diff, gene %in% rbps)
subset(dat.diff, gene %in% tfs.lst)
