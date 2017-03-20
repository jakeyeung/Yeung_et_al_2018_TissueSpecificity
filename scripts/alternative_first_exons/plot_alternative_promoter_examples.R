# 2016-12-07
# Jake Yeung
# plot_alternative_promoter_examples.R 

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")
library(ggrepel)
library(dplyr)
library(ggplot2)
library(hash)
library(reshape2)

source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/ModelStrToModel.R")

eps <- 1  # for log2 transform

# Load --------------------------------------------------------------------

do.filter <- TRUE
dataset <- "liverWTKO"
dataset <- "hogenesch"

tiss.filt <- c("Liver", "Kidney", "Cere", "BS", "Hypo")
if (dataset == "hogenesch"){
  load("Robjs/tpm.afe.avg.binary.Robj", verbose=T)
  load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)
  load("Robjs/tpm.merged.Robj", verbose=T); dat.bytranscript <- tpm.merged; rm(tpm.merged)
  load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj", v=T); fits.long.filt <- fits.best; rm(fits.best)
  tpm.afe.avg <- dplyr::rename(tpm.afe.avg, gene = gene_name)
  tpm.afe.avg <- dplyr::rename(tpm.afe.avg, transcript = transcript_id)
  tpm.afe.avg <- dplyr::rename(tpm.afe.avg, tpm = tpm.avg)
  
  tpm.counts <- subset(tpm.afe.avg, tissue == "Liver") %>%
    group_by(tissue, gene) %>%
    summarise(counts = length(transcript))
  tpm.counts <- tpm.counts[order(tpm.counts$counts, decreasing = T), ]
  
  counts.dic <- hash(as.character(tpm.counts$gene), tpm.counts$counts)
  
  if (do.filter){
    # tiss.filt <- c("Liver", "Kidney")
    tpm.afe.avg <- subset(tpm.afe.avg, tissue %in% tiss.filt)
    genes.filt <- as.character(tpm.afe.avg$gene)
    fits.long.filt <- subset(fits.long.filt, gene %in% genes.filt)
    dat.bytranscript <- subset(dat.bytranscript, tissue %in% tiss.filt)
  }
  dat.bytranscript <- dplyr::rename(dat.bytranscript, "gene" = gene_name)
  dat.bytranscript <- dplyr::rename(dat.bytranscript, "transcript" = transcript_id)
  dat.bytranscript$geno <- "WT"
  tstart <- 22
  tend <- 64
} else if (dataset == "liverWTKO"){
  load("Robjs/liver_kidney_atger_nestle/tpm.afe.avg.binary.Robj", v=T)
  load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
  load("Robjs/liver_kidney_atger_nestle/dat.bytranscript.bugfixed.Robj", v=T)
  load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
  dat.bytranscript <- StaggeredTimepointsLivKid(dat.bytranscript)
  dat.bytranscript$geno <- as.character(dat.bytranscript$geno)
  dat.bytranscript$geno <- gsub(pattern = "SV129", replacement = "WT", dat.bytranscript$geno)
  dat.bytranscript$geno <- gsub(pattern = "BmalKO", replacement = "BmalKO", dat.bytranscript$geno)
  dat.bytranscript$geno <- factor(as.character(dat.bytranscript$geno), levels = c("WT", "BmalKO"))
  dat.bytranscript.tg <- CollapseTissueGeno(dat.bytranscript)  # match fits.long.filt
  dat.long <- StaggeredTimepointsLivKid(dat.long)
  fits.long.filt <- subset(fits.long.filt, method == "g=1001")
  tstart <- 0
  tend <- 48
} else {
  stop("Dataset mustt be hogenesch or liverWTKO")
}


# Plot expression of examples ---------------------------------------------

jgenes <- c("Slc45a3", "Insig2", "Ddc", "Upp2", "Galnt11", "Gck", "Stat5b", "Plac8", "Ngef", "Psen2")
jtx <- c("ENSMUST00000027695", "ENSMUST00000177943")

dat.bytranscript$tissue <- factor(as.character(dat.bytranscript$tissue), levels = tiss.filt)

pdf(paste0("plots/alternative_exon_usage/gene_expression_across_transcripts.", dataset, ".pdf"))
jtx <- c("ENSMUST00000027695", "ENSMUST00000185387")
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == "Slc45a3" & tissue %in% c("Liver", "Kidney") & transcript %in% c("ENSMUST00000027695", "ENSMUST00000177943")), jtitle = "Slc45a3", log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == "Slc45a3" & tissue %in% c("Liver", "Kidney") & transcript %in% jtx), jtitle = "Slc45a3", log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))

# Ddc
jtx <- c("ENSMUST00000155690,ENSMUST00000134121,ENSMUST00000146359,ENSMUST00000109659",
         "ENSMUST00000066237,ENSMUST00000178704")
jgene <- "Ddc"
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & transcript %in% jtx), jtitle = jgene, log2.transform = FALSE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))
jtx <- c("ENSMUST00000155690,ENSMUST00000134121,ENSMUST00000146359,ENSMUST00000109659",
         "ENSMUST00000066237,ENSMUST00000178704")
jgene <- "Ddc"
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & transcript %in% jtx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))

# Insig2
jgene <- "Insig2"
jtx <- c("ENSMUST00000159085,ENSMUST00000159125,ENSMUST00000161818", 
         "ENSMUST00000186915,ENSMUST00000162582,ENSMUST00000160688,ENSMUST00000003818,ENSMUST00000160968")
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & transcript %in% jtx), jtitle = jgene, log2.transform = FALSE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))
jgene <- "Insig2"
jtx <- c("ENSMUST00000159085,ENSMUST00000159125,ENSMUST00000161818", 
         "ENSMUST00000186915,ENSMUST00000162582,ENSMUST00000160688,ENSMUST00000003818,ENSMUST00000160968")
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & transcript %in% jtx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))

# Gck
jtx <- c("ENSMUST00000102920", "ENSMUST00000125434")
jgene <- "Gck"
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & transcript %in% jtx), jtitle = jgene, log2.transform = FALSE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))
jtx <- c("ENSMUST00000102920", "ENSMUST00000125434")
jgene <- "Gck"
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & transcript %in% jtx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))

# Upp2
jtx <- c("ENSMUST00000135737,ENSMUST00000137007,ENSMUST00000071543", "ENSMUST00000059102,ENSMUST00000102755")
jgene <- "Upp2"
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & transcript %in% jtx), jtitle = jgene, log2.transform = FALSE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))
jtx <- c("ENSMUST00000135737,ENSMUST00000137007,ENSMUST00000071543", "ENSMUST00000059102,ENSMUST00000102755")
jgene <- "Upp2"
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & transcript %in% jtx), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))

for (jgene in jgenes){
  print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene), jtitle = jgene, log2.transform = FALSE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend)) 
  print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend)) 
}
# replot for WT only 
print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, geno == "WT" & gene == "Slc45a3" & transcript %in% c("ENSMUST00000027695", "ENSMUST00000177943")), jtitle = "Slc45a3", log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))
for (jgene in jgenes){
  print(PlotTpmAcrossTissuesWTKO(subset(dat.bytranscript, gene == jgene & geno == "WT"), jtitle = jgene, log2.transform = TRUE, transcript_id = "transcript", geno_id = "geno", tissue_id = "tissue", tstart = tstart, tend = tend))
}
dev.off()

