# 2016-07-09
# Make sitecounts table for N: use same parameters as in pLDA analysis

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

debug <- FALSE
# debug <- TRUE
print(paste("DEBUG:", debug))

tissue.spec.peaks <- TRUE

if (!debug){
  args <- commandArgs(trailingOnly=TRUE)
  distfilt <- as.numeric(args[1])
  jcutoff <- as.numeric(args[2])
  jcutoff.low <- 0.5
  jmethod <- args[3]
  jmodels <- c(args[4])  # "Liver_SV129"
} else {
  distfilt <- 40000
  jcutoff <- 3  # arbitrary
  jcutoff.low <- 0.5
  jmethod <- "g=1001"
  jmodels = c("Liver_SV129")
}

cleanup <- FALSE
writepeaks <- FALSE
# jmethod <- "BIC"

outdir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/liver_kidney_sitecounts_tissuespecpeaks_cutofflow"
dir.create(outdir, showWarnings = FALSE)
f <- paste0("sitecounts_enhancers.model_", jmodels, ".method_", jmethod, ".dist_", distfilt, ".cutoff_", jcutoff, ".cross_FALSE.", "tspeaks_", tissue.spec.peaks, ".cutofflow_", jcutoff.low, ".mat")
f.cross <- paste0("sitecounts_enhancers.model_", jmodels, ".method_", jmethod, ".dist_", distfilt, ".cutoff_", jcutoff, ".cross_TRUE.", "tspeaks_", tissue.spec.peaks, ".cutofflow_", jcutoff.low, ".mat")

outf.mat <- file.path(outdir, f)
outf.mat.cross <- file.path(outdir, f.cross)

if (file.exists(outf.mat) & file.exists(outf.mat.cross)){
  stop("File exists, stopping")
}

if (is.na(distfilt)) stop("Distfilt must be numeric")

print(paste("Distance filter:", distfilt))
print(paste("DHS signal zscore cutoff:", jcutoff))
print(paste("Method:", jmethod))
print(paste("Model:", jmodels))

saveplot <- FALSE
saverobj <- FALSE
amp.min <- 0

library(dplyr)
library(reshape2)
library(ggplot2)
library(hash)
library(penalizedLDA)

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/ShannonEntropy.R")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/PlotUCSC.R")  # CoordToBed
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotUCSC.R")
source("scripts/functions/LiverKidneyFunctions.R")



# Load --------------------------------------------------------------------


# Load --------------------------------------------------------------------

# save_N_on_posterior_cutoff_0.1.R saves Robj image. Here we laod it up
# Do Penalized LDA as before  ---------------------------------------------

# from multigene_analysis.play_with_parameters.R 
if (!exists("fits.best")){
  load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
  fits.best.hog <- fits.best
  
  load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.Robj", v=T)
  fits.best.orig <- fits.long.filt
  fits.best <- fits.long.filt; rm(fits.long.filt)
  fits.best <- subset(fits.best, method == jmethod)
  if (nrow(fits.best) == 0) stop("Method is wrong, fits.best is empty")
} 
if (!exists("dat.long")){
  load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
  dat.long.hog <- dat.long
  load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)
  dat.orig <- dat.long
  dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
  dat.long <- StaggeredTimepointsLivKid(dat.long)
} 
if (!exists("S.long")) load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
if (!exists("N.long.filt")){
  #   load("Robjs/N.long.livertwflat.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
  # load("Robjs/liver_kidney_atger_nestle/N.long.3wtmodules.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
  load("Robjs/liver_kidney_atger_nestle/N.long.all_genes.3wtmodules.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
}
# get fit from f24

load("Robjs/liver_kidney_atger_nestle/fits.bytiss.Robj", v=T)

fits.bytiss <- subset(fits.bytiss, tissue %in% jmodels & gene != "")
liver.amp <- hash(as.character(fits.bytiss$gene), fits.bytiss$amp)
liver.phase <- hash(as.character(fits.bytiss$gene), fits.bytiss$phase)


# Get genes and peaks -----------------------------------------------------


jgenes <- as.character(subset(fits.best, model %in% jmodels)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)


print(paste("Rhythmic genes:", length(jgenes)))
print(paste("Flat genes:", length(jgenes.flat)))


# print(paste("Rhythmic genes:", length(jgenes)))
# print(paste("Flat genes:", length(jgenes.flat)))

# get peaks near gene
S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
jpeaks.flat <- as.character(unique(S.sub.flat$peak))

print(paste("number of peaks in liver-specific rhythmic genes", length(jpeaks)))
print(paste("number of peaks in flat genes", length(jpeaks.flat)))

N.sub <- subset(N.long.filt, gene %in% jgenes & peak %in% jpeaks)  # should not have any missing peaks
N.sub.flat <- subset(N.long.filt, gene %in% jgenes.flat & peak %in% jpeaks.flat)

# Clean up ram ------------------------------------------------------------
if (cleanup){
  rm(S.long, N.long.filt)
}


# Filter peaks ------------------------------------------------------------

# tissue order as "Kidney", "Liver"

# tissue.spec.peaks = FALSE -> any Kidney and LIver above jcutoff
# tissue.spec.peaks = TRUE -> Liver above jcutoff, Kidney below jcutoff.low

if (!tissue.spec.peaks){
  print("Tissue agnostic peaks")
  S.sub.liverpeaks <- subset(S.sub, tissue %in% c("Kidney", "Liver")) %>%
    group_by(peak, gene) %>%  
    filter(max(zscore) > jcutoff)
} else {
  print("Tissue spec peaks")
  kidney.i <- 1
  liver.i <- 2
  S.sub.liverpeaks <- subset(S.sub, tissue %in% c("Kidney", "Liver")) %>%
    group_by(peak, gene) %>%  
    filter(zscore[liver.i] > jcutoff & zscore[kidney.i] < jcutoff.low)
}
mara.peaks <- unique(as.character(S.sub.liverpeaks$peak))

print(paste(length(mara.peaks), "mara peaks"))

# Create N matrix ---------------------------------------------------------

N.gene <- subset(N.sub, peak %in% mara.peaks) %>%
  group_by(gene, motif) %>%
  summarise(sitecount = sum(sitecount))

# check Mreg RORA
# check Tars RORA
# subset(N.gene, gene == "Tars" & motif == "RORA.p2")

mat.liver <- dcast(N.gene, formula = gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)

# write table
if (!file.exists(outf.mat) & debug == FALSE){
  print("Writing to:")
  print(outf.mat)
  write.table(mat.liver, file = outf.mat,
              quote = FALSE, sep = "\t", row.names = FALSE)
  
}



# Create N matrix cross product -------------------------------------------

rhyth.motifs <- GetTopMotifs("rhythmic")
rhyth.motifs <- c(rhyth.motifs, c("SRF.p3"))
rhyth.motifs <- rhyth.motifs[which( ! rhyth.motifs %in% c("ATF6.p2"))]
tissue.motifs <- GetTopMotifs("tissue")
tissue.motifs <- c(tissue.motifs, c("ATF5_CREB3.p2", "ATF6.p2"))
tissue.motifs <- tissue.motifs[which( ! tissue.motifs %in% c("SRF.p3"))]

# remove tissue motifs in rhyth
rhyth.motifs <- rhyth.motifs[which(!rhyth.motifs %in% intersect(rhyth.motifs, tissue.motifs))]

mat.rhyth <- subset(mat.liver, select = intersect(rhyth.motifs, colnames(mat.liver)))
mat.tiss <- subset(mat.liver, select = intersect(tissue.motifs, colnames(mat.liver)))
mat.rhythtiss <- CrossProductTwoSets(mat.rhyth, mat.tiss)

mat.liver.cross <- cbind(mat.liver, mat.rhythtiss)

if (!file.exists(outf.mat.cross) & debug == FALSE){
  print("Writing to:")
  print(outf.mat.cross)
  write.table(mat.liver.cross, file = outf.mat.cross,
              quote = FALSE, sep = "\t", row.names = FALSE)
  
}
