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
tiss <- c("Liver_SV129", "Kidney_SV129", "Liver_BmalKO", "Kidney_BmalKO")
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

# do only on subset of models??
jgenes <- unique(as.character(subset(fits.long.filt, model != "")$gene))
# jmods <- c("Liver_SV129", "Kidney_SV129")
# jgenes <- unique(as.character(subset(fits.long.filt, model %in% jmods)$gene))

# BEGIN: ANALYZE ON EACH TIME
tpm.afe.bysamp <- dat.bytranscript %>%
  group_by(gene, tissue, time, geno) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1)) %>%
  group_by(gene, transcript, tissue, geno, time) %>%
  summarise(tpm_norm.avg = mean(tpm_normalized),
            tpm_norm.var = var(tpm_normalized),
            tpm.avg = mean(tpm),
            tpm.var = var(tpm))
tpm.afe.bysamp <- subset(tpm.afe.bysamp, gene %in% unique(fits.long.filt$gene))
tpm.afe.bysamp$nprom <- sapply(as.character(tpm.afe.bysamp$gene), function(jgene) counts.dic[[jgene]])
tpm.afe.bysamp$amp <-  mapply(function(jtiss, jgene) is.rhyth.dic[[paste(jtiss, jgene, sep = ",")]], as.character(tpm.afe.bysamp$tissue), as.character(tpm.afe.bysamp$gene))

# Rename to just tissue, ignoring genotype
tpm.afe.bysamp$tissue.merge <- sapply(as.character(tpm.afe.bysamp$tissue), function(tiss) strsplit(tiss, "_")[[1]][[1]])

# amps of same tissue should be same (take maximum, if tissue in any genotype is rhythmic, then label it as rhythmic)
tpm.afe.bysamp <- tpm.afe.bysamp %>%
  group_by(gene, transcript, tissue.merge) %>%
  mutate(amp = max(amp))

# make tissue indivudal factors by adding geno and time
tpm.afe.bysamp$tissue <- paste(as.character(tpm.afe.bysamp$tissue), tpm.afe.bysamp$time, sep = "-")

tpm.gauss <- subset(tpm.afe.bysamp, nprom > 1 & gene %in% jgenes) %>%
  group_by(gene) %>%
  do(sigs = CalculateGaussianCenters(., transcript_id = "transcript"))

tpm.gauss2 <- subset(tpm.gauss, !is.na(sigs)) %>%
  group_by(gene) %>%
  do(CalculateGaussianDists(.))

# END: ANALYZE ON EACH TIME

genelist <- tpm.gauss2$gene
tpm.gauss <- cbind(as.data.frame(tpm.gauss2), as.data.frame(subset(tpm.gauss, gene %in% genelist, select = -gene)))

tpm.gauss <- tpm.gauss[order(tpm.gauss$center.dists, decreasing = TRUE), ]

head(data.frame(subset(tpm.gauss, select=-sigs)), n = 50)



# Testing single gene, commented out
# jgene <- "Mreg"
# # tpm.afe.bysamp <- subset(tpm.afe.bysamp, gene %in% unique(fits.long.filt$gene) & geno == "SV129")
# jtest <- subset(tpm.afe.bysamp, gene == jgene)
# jtest$tissue <- paste(jtest$tissue, jtest$time, sep = "-")
# out <- GetPromoterUsage(jtest, transcript_id = "transcript")
# 
# proms <- out$dat.mat.trans
# jtitle <- jgene
# plot(proms[, 2], proms[, 3], main = jtitle, xlab = "Promoter usage (1st component)", ylab = "Promoter usage (2nd component)")
# par(cex.axis=1, cex.lab=2, cex.main=2, cex.sub=1)
# 
# jcols <- vector(length = nrow(proms))
# jcols[which(proms$amp == 1)] <- "blue"
# jcols[which(proms$amp == 0)] <- "red"
# text(proms[, 2], proms[, 3], labels = proms$tissue, col = jcols)
# draw ellipse
# if (draw.ellipse){
#   try(lines(ellipse(tpm.gauss.sigs$sig1$cov, level = 0.5, centre = tpm.gauss.sigs$sig1$center), type='l', col = "blue"), silent = T)
#   try(lines(ellipse(tpm.gauss.sigs$sig2$cov, level = 0.5, centre = tpm.gauss.sigs$sig2$center), type='l', col = "red"), silent = T)
# }

# SANITY ------------------------------------------------------------------

# filter for similarly expressed genes (both high)

dat.min <- subset(dat.long, geno == "SV129") %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs)) %>%
  group_by(gene) %>%
  summarise(exprs.min = min(exprs.mean))

jcutoff <- 2
jgenes.min <- as.character(subset(dat.min, exprs.min >= jcutoff)$gene)

head(data.frame(subset(tpm.gauss, gene %in% jgenes.min, select=-sigs)), n = 50)

# TESTING
jgene <- "Slc45a3"
jgene <- "Gm14327"
jgene <- "Pnrc2"
PromoterSpacePlots(tpm.afe.bysamp, jgene = jgene, transcript_id = "transcript", use.weights = FALSE)
PlotGeneTissuesWTKO(subset(dat.long, gene == jgene), jtitle = jgene)

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


pdf("plots/alternative_exon_usage/liver_kidney.atger_nestle.bugfixed.bytime.pdf")
for (jgene in genes){
  print(jgene)
  jmodel <- as.character(subset(fits.long.filt, gene == jgene)$model[[1]])
  tx <- GetPromoterUsage(subset(tpm.afe.bysamp, gene == jgene), transcript_id = "transcript", get.entropy=FALSE, return.transcripts=TRUE)[[1]]$transcript
  print(PlotGeneTissuesWTKO(subset(dat.long, gene == jgene)) + ggtitle(paste(jgene, jmodel)))
  PromoterSpacePlots(tpm.afe.bysamp, jgene = jgene, transcript_id = "transcript", use.weights = FALSE)
  print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene & transcript %in% tx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
  print(PlotTpmAcrossTissues(subset(dat.bytranscript, gene == jgene), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript"))
}
dev.off()

# PromoterSpacePlots.nostics(subset(tpm.gauss, gene == jgene)$sigs[[1]], jgene, draw.ellipse = F)
# subset(fits.orig, gene == jgene)
