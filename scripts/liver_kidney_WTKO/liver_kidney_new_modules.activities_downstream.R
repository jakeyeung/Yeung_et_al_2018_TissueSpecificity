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
source("scripts/functions/NcondsFunctions.R")

# Load --------------------------------------------------------------------

plotdir <- "plots/mara_liver_kidney_modules"


load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == "g=1001")
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)
fits.long.filt$amp.avg <- mapply(GetAvgAmpFromParams, fits.long.filt$param.list, fits.long.filt$model)

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)
dat.orig <- dat.long
dat.long <- StaggeredTimepointsLivKid(dat.long)

# jmod <- "all"
jmod <- "Kidney_SV129,Kidney_BmalKO"
jmod <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO"
jmod <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129,Liver_BmalKO-Kidney_SV129,Liver_BmalKO,Kidney_BmalKO"
jmod <- "Liver_BmalKO"
jmod <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO"
jmod <- "Liver_SV129,Kidney_SV129"
jmod <- "many_modules_minrhyth.3"
jmod <- "Kidney_SV129"
jmod <- "Kidney_SV129,Kidney_BmalKO"
jmod <- "Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129"
jmod <- "many_modules_minrhyth.4"
jmod <- "many_modules_minrhyth.4.exclude_clockdriven_model"
jmod <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Liver_BmalKO.Kidney_SV129,Kidney_BmalKO"
outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
# outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.exclude_clockdriven_model"
outmain <- file.path(outbase, paste0("promoters.", jmod, ".g=1001"))
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.Gm129/promoters.", jmod, ".g=1001")
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod, ".g=1001")
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.removehic1/promoters.", jmod, ".g=1001")
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney./promoters.", jmod, ".g=1001")
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

# plot motifs
max.labs <- 15
jtitle <- ""
comp <- 1
eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = max.labs, jtitle = jtitle, peak.to.trough = TRUE, label.gene = c("bHLH_family.p2", "RORA.p2", "SRF.p3", "HSF1.2.p2", "TFAP2B.p2"))
print(eigens.act$u.plot + ggtitle(jmod))
multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)

# plot gene
if (grepl("^many_modules", jmod)){
  jn.rhyth <- as.numeric(strsplit(jmod, "\\.")[[1]][[2]])
  jmodels <- unique(as.character(subset(fits.long.filt, n.rhyth >= jn.rhyth)$model))
  if (jn.rhyth == 4){
    # exclude this model which is actually a clock-driven module
    excl.mod <- "Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO"
    jmodels <- jmodels[! jmodels %in% excl.mod]
  }
  # jmodels <- jmodels[!jmodels %in% "Liver_SV129,Liver_BmalKO"]
} else {
  jmodels <- jmod
}
genes.tw <- as.character(subset(fits.long.filt, model %in% jmodels)$gene)
s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = 1, label.n = 50, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
print(eigens$u.plot)

multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)


# Visualize top hits ------------------------------------------------------

tops <- head(unique(act.complex[order(Mod(act.complex$exprs.transformed), decreasing = TRUE), ]$gene), n = 35)

pdf(file.path(plotdir, paste0("top_hits_motifs.", jmod, ".pdf")))
for (jmotif in tops){
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
jmotifs <- head(as.character(act.tissdiff$gene), n = 70)

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
# Tcfap4
# Tcfe3

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
