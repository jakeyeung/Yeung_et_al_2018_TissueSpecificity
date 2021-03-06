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

# Functions ---------------------------------------------------------------

omega <- 2 * pi / 24
hr.shift <- 3
mrna.hl <- log(2) / (omega / tan(omega * hr.shift))  # convert hour shift to half-life 
mrna.hl.str <- signif(mrna.hl, digits = 2)

include.ko.prot <- FALSE


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

jmod <- "Kidney_SV129"
jmod.str <- gsub(",", "-", jmod)

gene.dat <- "Tfcp2"
gene.act <- "TFCP2"
gene.prot <- ""
jtiss <- "Kidney"

outbase <- "/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney"
outmain <- file.path(outbase, paste0("promoters.", jmod, ".g=1001"))
indir <- file.path(outmain, "atger_with_kidney.bugfixed")
source("scripts/functions/LoadActivitiesLong.R")
act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)

s.df <- GetAmpPhaseFromActivities(act.long, mrna.hl, jtiss = "Liver", jgeno = "WT")

act.shift <- act.long
act.shift$time <- act.long$time - hr.shift
act.shift$time <- sapply(act.shift$time, function(x) ifelse(x < 0, x + 48, x))


jwt.prot <- "CRY1"
pdf(paste0("/home/yeung/projects/tissue-specificity/plots/nuclear_proteomics_plots/mrna_activity_nuc_abund.IncludeKO.", include.ko.prot, ".Model.", jmod.str, ".NucProteinGeno.", jwt.prot, ".SingleDay.PhaseShift.", hr.shift, ".pdf"))
m <- PlotmRNAActivityProtein(dat.long, act.shift, prot.long = NA, gene.dat = gene.dat, gene.act = gene.act, gene.prot = gene.prot, jtiss = jtiss, single.day = FALSE) + theme(strip.text = element_blank())
m.day <- PlotmRNAActivityProtein(dat.long, act.shift, prot.long = NA, gene.dat = gene.dat, gene.act = gene.act, gene.prot = gene.prot, jtiss = jtiss, single.day = TRUE) + theme(strip.text = element_blank())
m.day.sine <- PlotmRNAActivityProtein(dat.long, s.df, prot.long = NA, gene.dat = gene.dat, gene.act = gene.act, gene.prot = gene.prot, jtiss = jtiss, single.day = TRUE, act.in.sine.cos = TRUE) + theme(strip.text = element_blank())
print(m)
print(m.day)
print(m.day.sine)
dev.off()