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
source("scripts/functions/ModelStrToModel.R")
source("/home/yeung/projects/posttranscriptional_regulation/functions/CosSineFunctions.R")

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
jmod <- "Kidney_SV129,Kidney_BmalKO"
jmod <- "many_modules_minrhyth.4"
jmod <- "many_modules_minrhyth.4.exclude_clockdriven_model"
jmod <- "Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129"
jmod <- "Kidney_SV129"
jmod <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Liver_BmalKO.Kidney_SV129,Kidney_BmalKO"
jtiss <- c("Liver_SV129", "Kidney_SV129", "Liver_BmalKO", "Kidney_BmalKO")


jmod <- "Liver_SV129"
jtiss <- c("Liver_SV129")

jmod <- "Liver_SV129,Liver_BmalKO"
jtiss <- c("Liver_SV129", "Liver_BmalKO")

outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
# outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.exclude_clockdriven_model"
outmain <- file.path(outbase, paste0("promoters.", jmod, ".g=1001"))
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.Gm129/promoters.", jmod, ".g=1001")
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod, ".g=1001")
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney.removehic1/promoters.", jmod, ".g=1001")
# outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney./promoters.", jmod, ".g=1001")
indir <- file.path(outmain, "atger_with_kidney.bugfixed")
source("scripts/functions/LoadActivitiesLong.R")
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)


# Plot module -------------------------------------------------------------

omega <- 2 * pi / 24

n <- 4

act.complex <- ProjectWithZscore(act.long, omega, n = 4)
# zscore.min <- 1
zscore.min <- 1.25
rhyth.tiss <- strsplit(jmod, ",")[[1]]
sig.motifs <- unique(as.character(subset(act.complex, zscore > zscore.min)$gene))

# act.complex <- act.long %>%
#   group_by(gene, tissue) %>%
#   do(ProjectWithZscore(., omega = omega, add.tissue = TRUE, propagate.errors = TRUE))
# head(print(data.frame(act.complex[order(act.complex$zscore, decreasing = TRUE), ])), n = 100)


# # Find z-score cutoff
# x <- seq(0, 3, length.out = 1000)
# y <- pnorm(x) - pnorm(-x)
# plot(x, y, type = "l")
# abline(h = 0.8)
# abline(v = 1.25)

rhyth.tiss <- strsplit(jmod, ",")[[1]]
# Filter out motifs whose amplitudes are smaller than amp.se 
sig.motifs <- unique(as.character(subset(act.complex, zscore > zscore.min)$gene))

print(sig.motifs)

s.act <- SvdOnComplex(subset(act.complex, gene %in% sig.motifs), value.var = "exprs.transformed")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
# jtitle <- gsub(pattern = "\\.", replacement = "\n", basename(indirmain))

# plot motifs
max.labs <- 25
jtitle <- ""
comp <- 1
eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = max.labs, jtitle = jtitle, peak.to.trough = TRUE, label.gene = c("bHLH_family.p2", "RORA.p2", "SRF.p3", "HSF1.2.p2", "TFAP2B.p2"))
print(eigens.act$u.plot + ggtitle(jmod))
multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)



# See if any corresponding proteins have rhythms --------------------------

source("scripts/functions/ProteomicsFunctions.R")
source("scripts/functions/GetTFs.R")

tfs <- GetTFs(split.commas = TRUE, get.motifs = TRUE, get.mat.only = TRUE)
prot.long <- LoadProteomicsData()
phospho.long <- LoadPhosphoData()

pdf("plots/nuclear_proteomics_plots/liver_WT_KO_hits.pdf")
for (jmotif in sig.motifs){
  print(jmotif)
  jgenes <- GetGenesFromMotifs(jmotif, tfs)
  for (jgene in jgenes){
    print(jgene)
    jsub <- subset(prot.long, gene == jgene)
    jsub.phospho <- subset(phospho.long, gene == jgene)
    # track phospho only (optional)
    if (nrow(jsub.phospho) == 0) next
      m <- PlotProteomics(jsub, jtitle = paste0("Nuc:", jgene))
      m.phos <- PlotProteomics(jsub.phospho, jtitle = paste0("Phos:", jgene))
      print(m)
      print(m.phos)
      for (jwt.prot in c("Cry", "Bmal")){
        nuc <- PlotmRNAActivityProtein(dat.long, act.long, prot.long, jgene, jmotif, jgene, jtiss = "Liver", dotsize = 3, themesize = 22, wt.prot = jwt.prot)
        phos <- PlotmRNAActivityProtein(dat.long, act.long, phospho.long, jgene, jmotif, jgene, jtiss = "Liver", dotsize = 3, themesize = 22, wt.prot = jwt.prot)
        print(nuc + ggtitle(paste0("Nuc:", jgene)))
        print(phos + ggtitle(paste0("Phos:", jgene)))
    }
  }
}
dev.off()


jgene <- "HOXA9_MEIS1.p2"
jgene <- "HIC1.p2"
jgene <- "HLF.p2"
jgene <- "bHLH_family.p2"
jgene <- "HSF1.2.p2"
jgene <- "SRF.p3"
jgene <- "ONECUT1.2.p2"
jgene <- "EGR1..3.p2"
jgene <- "HIC1.p2"
jgene <- "RORA.p2"
PlotActivitiesWithSE(subset(act.long, gene == jgene), jtitle = "", showSE = TRUE) + ggtitle(jgene)

# check if it is centered??
act.sum <- act.long %>%
  group_by(gene) %>%
  summarise(exprs = sum(exprs))

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
