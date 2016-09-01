# 2016-08-16
# Jake Yeung
# explore_proteomics_data.R

rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)

source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/ProteomicsFunctions.R")

# Functions ---------------------------------------------------------------


prot.long <- LoadProteomicsData()


# Get activities ----------------------------------------------------------


load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == "g=1001")
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)
fits.long.filt$amp.avg <- mapply(GetAvgAmpFromParams, fits.long.filt$param.list, fits.long.filt$model)

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)
dat.orig <- dat.long
dat.long <- StaggeredTimepointsLivKid(dat.long)

jmod <- "Liver_SV129,Liver_BmalKO"

outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
outmain <- file.path(outbase, paste0("promoters.", jmod, ".g=1001"))
indir <- file.path(outmain, "atger_with_kidney.bugfixed")
source("scripts/functions/LoadActivitiesLong.R")
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = FALSE)


# Plot proteomics, mRNA, and activity all on same plot --------------------

gene.lst <- list()
gene.dat <- "Stat2"
gene.act <- "STAT2.4.6.p2"
gene.prot <- "Stat2"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

gene.dat <- "Elf2"
gene.act <- "ELF1.2.4.p2"
gene.prot <- "Elf2"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

gene.dat <- "Junb"
gene.act <- "JUN.p2"
gene.prot <- "Junb"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

gene.dat <- "Creb3"
gene.act <- "ATF5_CREB3.p2"
gene.prot <- "Crebbp"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

gene.dat <- "Egr1"
gene.act <- "EGR1..3.p2"
gene.prot <- "Egr1"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

gene.dat <- "Mafb"
gene.act <- "MAFB.p2"
gene.prot <- "Mafb"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

gene.dat <- "Tfdp1"
gene.act <- "TFDP1.p2"
gene.prot <- "Tfdp1"

# gene.dat <- "Nrf1"
# gene.act <- "NRF1.p2"
# gene.prot <- "Nrf1"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

# setEPS()
# postscript("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.eps")
# postscript("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.ps")

pdf("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.pdf")
for (genes in gene.lst){
  gene.dat <- genes[1]; gene.act <- genes[2]; gene.prot <- genes[3]
  jtest <- PlotmRNAActivityProtein(dat.long, act.long, prot.long, gene.dat, gene.act, gene.prot, jtiss = "Liver", dotsize = 3, themesize = 22)
  jtest.both <- PlotmRNAActivityProtein(dat.long, act.long, prot.long, gene.dat, gene.act, gene.prot, jtiss = "both", dotsize = 2, themesize = 14)
  print(jtest + ggtitle(gene.dat))
  print(jtest.both + ggtitle(gene.dat) + theme(strip.text = element_blank()))
}
dev.off()

# Explore rhythms in proteomics data --------------------------------------

jgene <- "Stat2"
jgene <- "Srfbp1"
jgene <- "Nfil3"
jgene <- "Nr3c1"
jgene <- "Nr3c2"
jgene <- "Atf2"
jgene <- "Tfdp1"
jgene <- "Hmga2"
PlotProteomics(subset(prot.long, gene == jgene)) + ggtitle(jgene)

jgenes <- c("Stat2", "Mafb", "Junb", "E2f5", "Creb3", "Elf2", "Tgif1", "Egr1", "Pou2f1", "Hes6", "Tfdp1", "Crebbp")

jgene <- "Onecut1"
jgene <- "Onecut2"
jgene <- "Nfya"
jgene <- "Nfyc"
jgene <- "Hic1"
# print(PlotProteomics(jsub, jtitle = jgene))

for (jgene in jgenes){
  print(jgene)
  jsub <- subset(prot.long, gene == jgene)
  if (nrow(jsub) > 0 & !all(is.na(jsub$rel.abund))){
    print(PlotProteomics(jsub, jtitle = jgene))
  } else {
    print(paste("Skipping", jgene))
  }
}
