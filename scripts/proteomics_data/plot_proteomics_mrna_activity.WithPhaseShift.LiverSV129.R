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
source("scripts/functions/HalfLifeFunctions.R")
source("scripts/functions/CosSineFunctions.R")


# Set constants -----------------------------------------------------------

omega <- 2 * pi / 24
hr.shift <- 3
mrna.hl <- log(2) / (omega / tan(omega * hr.shift))  # convert hour shift to half-life 
mrna.hl.str <- signif(mrna.hl, digits = 2)

include.ko.prot <- FALSE


# Functions ---------------------------------------------------------------


prot.long <- LoadProteomicsData()

if (!include.ko.prot){
  prot.long <- subset(prot.long, geno == "WT")
}


# Get activities ----------------------------------------------------------


load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == "g=1001")
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)
fits.long.filt$amp.avg <- mapply(GetAvgAmpFromParams, fits.long.filt$param.list, fits.long.filt$model)

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)
dat.orig <- dat.long
dat.long <- StaggeredTimepointsLivKid(dat.long)

jmod <- "Liver_SV129"
jmod.str <- gsub(",", "-", jmod)

outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
outmain <- file.path(outbase, paste0("promoters.", jmod, ".g=1001"))
indir <- file.path(outmain, "atger_with_kidney.bugfixed")
source("scripts/functions/LoadActivitiesLong.R")
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = FALSE)


# Plot proteomics, mRNA, and activity all on same plot --------------------

gene.lst <- list()

gene.dat <- "Elf1"
gene.act <- "ELF1.2.4.p2"
gene.prot <- "Elf1"
gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

s.df <- GetAmpPhaseFromActivities(act.long, mrna.hl, jtiss = "Liver", jgeno = "WT")

act.shift <- act.long
act.shift$time <- act.long$time - hr.shift
act.shift$time <- sapply(act.shift$time, function(x) ifelse(x < 0, x + 48, x))

# setEPS()
# postscript("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.eps")
# postscript("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.ps")

# jwt.prot <- "Bmal" 
for (jwt.prot in c("Bmal", "Cry")){
  # pdf(paste0("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.Model.", jmod, ".NucProteinGeno.", jwt.prot, ".SingleDay.PhaseShift.", hr.shift, ".pdf"))
  pdf(paste0("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.IncludeKO.", include.ko.prot, ".Model.", jmod.str, ".NucProteinGeno.", jwt.prot, ".SingleDay.PhaseShift.", hr.shift, ".pdf"))
  for (genes in gene.lst){
    gene.dat <- genes[1]; gene.act <- genes[2]; gene.prot <- genes[3]
    # jtest <- PlotmRNAActivityProtein(dat.long, act.long, prot.long, gene.dat, gene.act, gene.prot, jtiss = "Liver", dotsize = 3, themesize = 22, wt.prot = jwt.prot)
    # jtest.both <- PlotmRNAActivityProtein(dat.long, act.long, prot.long, gene.dat, gene.act, gene.prot, jtiss = "both", dotsize = 2, themesize = 14, wt.prot = jwt.prot)
    jtest.singleday.ampphase <- PlotmRNAActivityProtein(dat.long, s.df, prot.long = prot.long, gene.dat = gene.dat,  gene.act = gene.act, gene.prot = gene.prot, jtiss = "Liver", dotsize = 3, themesize = 22, wt.prot = jwt.prot, line.for.protein = TRUE, act.in.sine.cos = TRUE, single.day = TRUE) + theme(strip.text = element_blank())
    jtest.singleday <- PlotmRNAActivityProtein(dat.long, act.shift, prot.long = prot.long, gene.dat = gene.dat,  gene.act = gene.act, gene.prot = gene.prot, jtiss = "Liver", dotsize = 3, themesize = 22, wt.prot = jwt.prot, line.for.protein = TRUE, act.in.sine.cos = FALSE, single.day = TRUE) + theme(strip.text = element_blank())
    # print(jtest + ggtitle(gene.dat))
    # print(jtest.both + ggtitle(gene.dat) + theme(strip.text = element_blank()))
    print(jtest.singleday + ggtitle(gene.dat))
    print(jtest.singleday.ampphase + ggtitle(gene.dat))
  }
  dev.off()
}
  
  

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
