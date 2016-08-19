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
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
dat.orig <- dat.long
dat.long <- StaggeredTimepointsLivKid(dat.long)

# jmod <- "all"
jmod <- "Kidney_SV129,Kidney_BmalKO"
jmod <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO"
jmod <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129,Liver_BmalKO-Kidney_SV129,Liver_BmalKO,Kidney_BmalKO"
jmod <- "Liver_BmalKO"
jmod <- "Liver_SV129,Kidney_SV129"
jmod <- "many_modules_minrhyth.4"
jmod <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO"
# outmain <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.Kidney_SV129,Kidney_BmalKO.g=1001"
outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod, ".g=1001")
indir <- file.path(outmain, "atger_with_kidney.bugfixed")
source("scripts/functions/LoadActivitiesLong.R")
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = FALSE)


# Plot module -------------------------------------------------------------

omega <- 2 * pi / 24
act.complex <- act.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))

s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")


jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# jtitle <- gsub(pattern = "\\.", replacement = "\n", basename(indirmain))

max.labs <- 30
jtitle <- ""
comp <- 1
eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = max.labs, jtitle = jtitle, peak.to.trough = TRUE, label.gene = c("bHLH_family.p2", "RORA.p2", "SRF.p3", "HSF1.2.p2"))
print(eigens.act$u.plot + ggtitle(jmod))
# multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)

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

act.tissdiff <- act.tissmean %>%  group_by(gene) %>%
  summarise(exprs.diff = diff(exprs)) %>%
  # arrange(desc(abs(exprs.diff)))
  arrange(exprs.diff)

print(head(act.tissdiff))  # MYBL2.p2 etc
jmotifs <- head(as.character(act.tissdiff$gene), n = 35)
plotdir <- "plots/mara_liver_kidney_modules"
pdf(file.path(plotdir, paste0("tissue_spec_motifs.", jmod, ".pdf")))
# jmotif <- "HSF1.2.p2"
for (jmotif in jmotifs){
  print(PlotActivitiesWithSE.rnaseq(subset(act.long, gene == jmotif)))
  tfs <- GetTFs(get.mat.only = TRUE)
  for (g in GetGenesFromMotifs(jmotif, tfs)){
    print(g)
    jsub <- subset(dat.orig, gene == g)
    if (nrow(jsub) > 0){
      print(PlotGeneTissuesWTKO(jsub, jtitle = g))
    } else {
      print(paste("Could not find gene:", g))
    }
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
