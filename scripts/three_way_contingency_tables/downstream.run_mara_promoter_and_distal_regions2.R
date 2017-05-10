# 2017-05-03
# Jake Yeung
# Downstream MARA: can we identify promoter-specific and distal-specific motifs?

rm(list=ls())

library(ggplot2)
library(dplyr)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/CosSineFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/GetTFs.R")

FixTissueGenoNames <- function(act.s){
  act.s$tissue <- sapply(as.character(act.s$tissue), function(x) strsplit(x, "_")[[1]][[1]])
  # act.s$tissue <- factor(as.character(act.s$tissue), levels = c("Liver_SV129", "Liver_BmalKO", "Kidney_SV129", "Kidney_BmalKO"))
  act.s$tissue <- factor(as.character(act.s$tissue), levels = c("Liver", "Kidney"))
  act.s$geno <- gsub("SV129", "WT", act.s$geno)
  act.s$geno <- gsub("BmalKO", "Bmal1 KO", act.s$geno)
  act.s$geno <- factor(as.character(act.s$geno), levels = c("WT", "Bmal1 KO"))
  return(act.s)
}

# Load tables -------------------------------------------------------------

hr.shift <- 3
zscore.min <- 1.25
n <- 4
omega <- 2 * pi / 24
jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)

# Load --------------------------------------------------------------------

# promoter motifs

jmod <- "Liver_SV129,Liver_BmalKO"
jmod <- "Liver_SV129"
jmodstr <- gsub(",", "-", jmod)
jtiss <- strsplit(jmod, "_")[[1]][[1]]
print(paste("Finding regulators for:", jmod))
genes.mod <- unique(as.character(subset(fits.long.filt, model == jmod)$gene))

maraoutdir <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod, ".g=1001/atger_with_kidney.bugfixed")
act.s <- LoadActivitiesLongKidneyLiver(maraoutdir, shorten.motif.name = TRUE)
act.s.orig <- act.s  # for ProjectingWithZscore

act.s <- FixTissueGenoNames(act.s)
# act.s$tissue <- sapply(as.character(act.s$tissue), function(x) strsplit(x, "_")[[1]][[1]])
# # act.s$tissue <- factor(as.character(act.s$tissue), levels = c("Liver_SV129", "Liver_BmalKO", "Kidney_SV129", "Kidney_BmalKO"))
# act.s$tissue <- factor(as.character(act.s$tissue), levels = c("Liver", "Kidney"))
# act.s$geno <- gsub("SV129", "WT", act.s$geno)
# act.s$geno <- gsub("BmalKO", "Bmal1 KO", act.s$geno)
# act.s$geno <- factor(as.character(act.s$geno), levels = c("WT", "Bmal1 KO"))

act.s.complex <- ProjectWithZscore(act.s.orig, omega, n)
sig.motifs <- unique(as.character(subset(act.s.complex, zscore > zscore.min)$gene))
s.act <- SvdOnComplex(subset(act.s.complex, gene %in% sig.motifs), value.var = "exprs.transformed")
jmod.label <- gsub(pattern = ",", replacement = "-", jmod)
# just shift by 3 hours
act.s.shift <- act.s
act.s.shift$time <- act.s$time - hr.shift
act.s.shift$time <- sapply(act.s.shift$time, function(x) ifelse(x < 0, x + 48, x))

# plot example

PlotActivitiesWithSE.wtko(subset(act.s.shift, gene == "RORA")) + theme(aspect.ratio = 1, strip.text = element_text())
PlotActivitiesWithSE.wtko(subset(act.s.shift, gene == "AR")) + theme(aspect.ratio = 1, strip.text = element_text())
PlotActivitiesWithSE.wtko(subset(act.s.shift, gene == "bHLH_family")) + theme(aspect.ratio = 1, strip.text = element_text())
PlotActivitiesWithSE.wtko(subset(act.s.shift, gene == "bHLH_family")) + theme(aspect.ratio = 1, strip.text = element_text())


# distal motifs
maraoutdir.distal <- "/home/yeung/data/tissue_specificity/mara_results/mara_outputs.weight.0.sql.TRUE.mod.Liver_SV129.dhscutoff.3.0/center.TRUE.weight.0.sql.TRUE.mod.Liver_SV129.dhscutoff.3.0/centered.TRUE.mod.Liver_SV129.weightcutoff.0"
cutstr <- "2.5.0"
cutstr <- "2.0"
cutstr <- "3.0"
cutstr <- "1.5.0"
maraoutdir.distal <- paste0("/home/yeung/data/tissue_specificity/mara_results/mara_outputs.weight.0.sql.TRUE.mod.", jmodstr, ".dhscutoff.", cutstr, ".final/center.TRUE.weight.0.sql.TRUE.mod.", jmodstr, ".dhscutoff.", cutstr, ".final/centered.TRUE.mod.", jmodstr, ".weightcutoff.0")
act.s.distal <- LoadActivitiesLong(maraoutdir.distal, shorten.motif.name = TRUE, make.cnames = FALSE)
act.s.distal <- MakeCnamesLivKidWTKO(act.s.distal)
act.s.distal <- FixTissueGenoNames(act.s.distal)
# act.s.distal$tissue <- factor(as.character(act.s.distal$tissue), levels = c("Liver_SV129", "Liver_BmalKO", "Kidney_SV129", "Kidney_BmalKO"))
act.s.shift.distal <- act.s.distal
act.s.shift.distal$time <- act.s.distal$time - hr.shift
act.s.shift.distal$time <- sapply(act.s.shift.distal$time, function(x) ifelse(x < 0, x + 48, x))

# PlotActivitiesWithSE.wtko(subset(act.s.shift.distal, gene == "RORA")) + theme(aspect.ratio = 1, strip.text = element_text())
# PlotActivitiesWithSE.wtko(subset(act.s.shift.distal, gene == "AR")) + theme(aspect.ratio = 1, strip.text = element_text())
# PlotActivitiesWithSE.wtko(subset(act.s.shift.distal, gene == "bHLH_family")) + theme(aspect.ratio = 1, strip.text = element_text())

tiss.filt <- c("Liver_SV129", "Liver_BmalKO")
jmotif <- "ELF1.2.4"
jmotif <- "NR3C1"
jmotif <- "HOX.A4.D4."
jmotif <- "FOXA2"
jmotif <- "ELF1.2.4"
jmotif <- "MAFB"
jmotif <- "ESR1"
jmotif <- "EGR1..3"
jmotif <- "CEBPA.B_DDIT3"
jmotif <- "bHLH_family"
jmotif <- "RORA"

if (jmod == "Liver_SV129"){
  hits <- c("ELF1.2.4", "bHLH_family", "RORA", "AR", "JUN", "HAND1.2", "ARNT_ARNT2_BHLHB2_MAX_MYC_USF1")
} else if (jmod == "Liver_SV129,Liver_BmalKO"){
  hits <- c("MAFB", "EGR1..3", "CEBPA.B_DDIT3", "RXRG_dimer", "ESR1")
}

tfs <- GetTFs(get.mat.only = TRUE, shorten.motif.name = TRUE)
rownames(tfs) <- make.names(rownames(tfs))

pdf(paste0("/home/yeung/projects/tissue-specificity/plots/mara_liver_kidney_modules_on_liverDHS.summary/mod.", jmodstr, ".cutoff.", cutstr, ".promoter_distal_hits.singleday.pdf"), useDingbats = FALSE)
for (h in hits){
  for (jshowSE in c(TRUE, FALSE)){
    m.p <- PlotActivitiesWithSE.wtko(subset(act.s.shift, gene == h), showSE = jshowSE, single.day=TRUE) + ggtitle(paste(h, "PromoterOnly"))  + theme(aspect.ratio = 1, strip.text = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    m.d <- PlotActivitiesWithSE.wtko(subset(act.s.shift.distal, gene == h), showSE = jshowSE, single.day=TRUE) + ggtitle(paste(h, "Distal", cutstr)) + theme(aspect.ratio = 1, strip.text = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m.p)
    print(m.d)
    multiplot(m.p + theme_bw(10) + theme(aspect.ratio = 1, strip.text = element_blank()), m.d + theme_bw(10) + theme(aspect.ratio = 1, strip.text = element_blank()), cols = 1)
  }
  # PlotActivitiesWithSE.wtko(subset(act.s.shift, gene == "RORA"), showSE = TRUE) + theme(aspect.ratio = 1, strip.text = element_text())
  # PlotActivitiesWithSE.wtko(subset(act.s.shift, gene == "RORA"), showSE = FALSE) + theme(aspect.ratio = 1, strip.text = element_text())
  genes.all <- unlist(sapply(h, GetGenesFromMotifs, tfs))
  for (g in genes.all){
    print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == g), jtitle = g, single.day = TRUE))
  }
}
dev.off()

