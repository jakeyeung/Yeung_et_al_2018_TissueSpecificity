# 2016-07-01
# Jake Yeung
# Alternative promoter usage

rm(list=ls())

library(dplyr)
library(ggplot2)
library(hash)
library(mvtnorm)
library(reshape2)

source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/NcondsFunctions.R")

eps <- 1  # for log2 transform

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.bytranscript.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)

fits.orig <- fits.long.filt
fits.long.filt <- subset(fits.long.filt, method == "g=1001")

# Annotate fits
fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)
fits.long.filt$amp.avg <- mapply(GetAvgAmpFromParams, fits.long.filt$param.list, fits.long.filt$model)
fits.long.filt$phase.sd <- mapply(GetSdPhaseFromParams, fits.long.filt$param.list, fits.long.filt$model)
fits.long.filt$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.long.filt$param.list, fits.long.filt$model)
fits.long.filt$phase.avg <- mapply(GetPhaseFromParams, fits.long.filt$param.list, fits.long.filt$model)


dat.bytranscript <- StaggeredTimepointsLivKid(dat.bytranscript)
dat.bytranscript <- CollapseTissueGeno(dat.bytranscript)  # match fits.long.filt


dat.long <- StaggeredTimepointsLivKid(dat.long)

tpm.afe.avg <- dat.bytranscript %>%
  group_by(gene, tissue, time, geno) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1)) %>%
  group_by(gene, transcript, tissue, geno) %>%
  summarise(tpm_norm.avg = mean(tpm_normalized),
            tpm_norm.var = var(tpm_normalized),
            tpm.avg = mean(tpm),
            tpm.var = var(tpm))


# Positive control --------------------------------------------------------

jgene <- "Ddc"
jgene <- "Slc45a3"

jsub <- subset(dat.bytranscript, gene == jgene & geno == "SV129")

PlotTpmAcrossTissues(jsub, jtitle = jgene, log2.transform = FALSE, transcript_id = "transcript")
PlotGeneTissuesWTKO(subset(dat.long, gene == jgene))


# Count promoters ---------------------------------------------------------
tpm.counts <- subset(tpm.afe.avg, tissue == "Liver_SV129") %>%
  group_by(tissue, gene) %>%
  summarise(counts = length(transcript))
tpm.counts <- tpm.counts[order(tpm.counts$counts, decreasing = T), ]

counts.dic <- hash(as.character(tpm.counts$gene), tpm.counts$counts)
tpm.afe.avg$nprom <- sapply(as.character(tpm.afe.avg$gene), function(jgene) counts.dic[[jgene]])

# hash it
tiss <- c("Liver_SV129", "Kidney_SV129")
keys.df <- fits.long.filt %>%
  group_by(gene) %>%
  do(GetGeneModelKeys(., tiss))
keys <- paste(keys.df$tissue, keys.df$gene, sep = ",")
vals <- as.numeric(keys.df$is.rhyth)
is.rhyth.dic <- hash(keys, vals)

# annotate it
tpm.afe.avg <- subset(tpm.afe.avg, gene %in% unique(fits.long.filt$gene) & geno == "SV129")
# 1 or 0 to signify whether it is rhythmic
tpm.afe.avg$amp <- mapply(function(jtiss, jgene) is.rhyth.dic[[paste(jtiss, jgene, sep = ",")]], as.character(tpm.afe.avg$tissue), as.character(tpm.afe.avg$gene))

# do only for Liver or Kidney rhythmic genes??
jmods <- c("Liver_SV129", "Kidney_SV129")
jgenes <- unique(as.character(subset(fits.long.filt, model %in% jmods)$gene))

tpm.gauss <- subset(tpm.afe.avg, nprom > 1 & gene %in% jgenes) %>%
  group_by(gene) %>%
  do(sigs = CalculateGaussianCenters(., transcript_id = "transcript"))

tpm.gauss2 <- subset(tpm.gauss, !is.na(sigs)) %>%
  group_by(gene) %>%
  do(CalculateGaussianDists(.))

genelist <- tpm.gauss2$gene
tpm.gauss <- cbind(as.data.frame(tpm.gauss2), as.data.frame(subset(tpm.gauss, gene %in% genelist, select = -gene)))

tpm.gauss <- tpm.gauss[order(tpm.gauss$center.dists, decreasing = TRUE), ]

head(data.frame(subset(tpm.gauss, select=-sigs)), n = 50)

# SANITY ------------------------------------------------------------------

# filter for similarly expressed genes (both high)

dat.min <- subset(dat.long, geno == "SV129") %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs)) %>%
  group_by(gene) %>%
  summarise(exprs.min = min(exprs.mean))

jcutoff <- 2
jgenes <- as.character(subset(dat.min, exprs.min >= 2)$gene)

head(data.frame(subset(tpm.gauss, gene %in% jgenes, select=-sigs)), n = 50)

jgene <- "Npnt"
jgene <- "Osgin1"
jgene <- "Zranb1"
jgene <- "Aox3"
jgene <- "Upp2"
jgene <- "9530068E07Rik"
jgene <- "Prkd3"
jgene <- "Sh3bgrl2"
jgene <- "Por"
jgene <- "Aqp9"
jgene <- "Npnt"
jgene <- "Wdtc1"
jgene <- "Zranb1"
jgene <- "Psen2"
jgene <- "Prkd3"
jgene <- "Insig2"
jgene <- "Upp2"
jgene <- "Ddc"
jgene <- "Usp2"
jgene <- "Slc45a3"

genes <- as.character(head(data.frame(subset(tpm.gauss, gene %in% jgenes, select=-sigs)), n = 50)$gene)

pdf("plots/alternative_exon_usage/liver_kidney.atger_nestle.bugfixed.pdf")
for (jgene in genes){
  print(jgene)
  tx <- GetPromoterUsage(subset(tpm.afe.avg, gene == jgene), transcript_id = "transcript", get.entropy=FALSE, return.transcripts=TRUE)[[1]]$transcript
  print(PlotGeneTissuesWTKO(subset(dat.long, gene == jgene)) + ggtitle(jgene))
  print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & geno == "SV129" & transcript %in% tx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
  print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & geno == "SV129"), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
}
dev.off()

# PromoterSpacePlots.nostics(subset(tpm.gauss, gene == jgene)$sigs[[1]], jgene, draw.ellipse = F)
# subset(fits.orig, gene == jgene)
