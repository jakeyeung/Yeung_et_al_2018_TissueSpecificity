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


# Load --------------------------------------------------------------------

inf <- "/home/shared/nuclear_proteomics/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.OneGenePerLine.txt"

prot <- read.table(inf, header = TRUE, sep = "\t")


# Make long ---------------------------------------------------------------

wt.sampnames <- paste("ZT", sprintf("%02d", seq(0, 45, 3)), ".WT", sep = "")
ko.sampnames <- paste("ZT", sprintf("%02d", seq(0, 18, 6)), ".Bmal.WT", sep = "")
fit.sampnames <- c("mean", "amp", "relamp", "phase", "pval", "qv", "amp.12h", "relamp.12h", "phase.12h", "pval.12h", "qv.12h")

prot.long.wt <- melt(prot, id.vars = "Gene.names", measure.vars = wt.sampnames, variable.name = "samp", value.name = "rel.abund")
prot.long.bmalko <- melt(prot, id.vars = "Gene.names", measure.vars = ko.sampnames, variable.name = "samp", value.name = "rel.abund")
fit.prot.wt <- subset(prot, select = c("Gene.names", fit.sampnames))
# fit.prot.wt <- melt(prot, id.vars = "Gene.names", measure.vars = fit.sampnames)

prot.long.wt$time <- GetTimeFromSamp(as.character(prot.long.wt$samp))
prot.long.wt$geno <- GetGenoFromSamp(as.character(prot.long.wt$samp))

prot.long.bmalko$time <- GetTimeFromSamp(as.character(prot.long.bmalko$samp))
prot.long.bmalko$geno <- GetGenoFromSamp(as.character(prot.long.bmalko$samp))

# merge
prot.long <- rbind(prot.long.wt, prot.long.bmalko)

# change Gene.names to gene
colnames(fit.prot.wt)[which(colnames(fit.prot.wt) == "Gene.names")] <- "gene"
colnames(prot.long)[which(colnames(prot.long) == "Gene.names")] <- "gene"


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

gene.dat <- "Tfdp1"
gene.act <- "TFDP1.p2"
gene.prot <- "Tfdp1"

# gene.dat <- "Nrf1"
# gene.act <- "NRF1.p2"
# gene.prot <- "Nrf1"

gene.lst[[gene.dat]] <- c(gene.dat, gene.act, gene.prot)

pdf("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.pdf")
for (genes in gene.lst){
  gene.dat <- genes[1]; gene.act <- genes[2]; gene.prot <- genes[3]
  jtest <- PlotmRNAActivityProtein(dat.long, act.long, prot.long, gene.dat, gene.act, gene.prot)
  print(jtest + ggtitle(gene.dat))
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
