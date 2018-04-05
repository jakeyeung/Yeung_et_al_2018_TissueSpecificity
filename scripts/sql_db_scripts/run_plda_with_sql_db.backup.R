# Jake Yeung
# run_plda_with_sql_db.backup.R
#  
# 2016-09-14

# 2016-04-21
# Jake Yeung
# pLDA was done before with cutoff at 0.5, but no real reason to do that. Redo with 0.1!

# 2016-04-21
# Jake Yeung
# pLDA was done before with cutoff at 0.5, but no real reason to do that. Redo with 0.1!

rm(list=ls())

start <- Sys.time()

setwd("/home/yeung/projects/tissue-specificity")

args <- commandArgs(trailingOnly=TRUE)
distfilt <- as.numeric(args[1])
jcutoff <- as.numeric(args[2])
jcutoff.low <- as.numeric(args[3])
distfilt <- 40000
jcutoff <- 3  # arbitrary
jcutoff.low <- 0  # arbitrary
# jcutoff <- 2  # arbitrary
# jcutoff <- 3  # arbitrary
cleanup <- FALSE
writepeaks <- FALSE
jmethod <- "g=1001"

if (is.na(distfilt)) stop("Distfilt must be numeric")

print(paste("Distance filter:", distfilt))
print(paste("DHS signal zscore cutoff:", jcutoff))

saveplot <- FALSE
outdir <- "plots/penalized_lda/liver_kidney_wtko.sql"
dir.create(outdir)
outf <- paste0(outdir, "2D.posterior.multigene.distfilt.kidneyWTKO.", distfilt, ".cutoff.", jcutoff, ".cutofflow", jcutoff.low, ".method.", jmethod, ".pdf")
amp.min <- 0

jmodels <- c("Kidney_SV129")
# jmodels <- c("Liver_SV129")
if (jmodels == "Kidney_SV129"){
  rhyth.tiss <- c("Kidney")
  flat.tiss <- c("Liver")
} else if (jmodels == "Liver_SV129"){
  rhyth.tiss <- c("Liver")
  flat.tiss <- c("Kidney")
}
outfile.robj <- paste0("Robjs/liver_kidney_atger_nestle/penalized_lda_mats.posterior.model.", jmodels[[1]], ".distfilt.", distfilt, ".", paste(rhyth.tiss, sep = "_"), ".cutoff", jcutoff, ".method.", jmethod, ".Robj")
if (file.exists(outfile.robj)) stop(paste0(outfile.robj, " exists. Exiting"))

library(dplyr)
library(reshape2)
library(ggplot2)
library(hash)
library(penalizedLDA)

source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")
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


# Functions ---------------------------------------------------------------

colMax <- function(dat){
  return(apply(dat, 1, max))
}

# Load --------------------------------------------------------------------

# Load up sitecounts from sql database
inf <- "/home/shared/sql_dbs/closestbed_multiple_genes.genomewide.merged.motifindexed.sqlite3"
motevo.tbl <- LoadDatabase(inf)

# from multigene_analysis.play_with_parameters.R 
if (!exists("fits.best")){
  load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
  fits.best <- fits.long.filt; rm(fits.long.filt)
  fits.best <- subset(fits.best, method == jmethod)
} 
if (!exists("dat.long")){
  load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
  dat.orig <- dat.long
  dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
  dat.long <- StaggeredTimepointsLivKid(dat.long)
} 
if (!exists("S.long")) load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
# load N.long.filt after we defined thhe genes we're looking for 



# Get genes and peaks -----------------------------------------------------


jgenes <- as.character(subset(fits.best, model %in% jmodels)$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)


print(paste("Rhythmic genes:", length(jgenes)))
print(paste("Flat genes:", length(jgenes.flat)))


# get peaks near gene
S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
jpeaks.flat <- as.character(unique(S.sub.flat$peak))

print(paste("number of peaks in liver-specific rhythmic genes", length(jpeaks)))
print(paste("number of peaks in flat genes", length(jpeaks.flat)))


N.long.filt.query <- filter(motevo.tbl, gene %in% jgenes)  # peaks are not indexed, so dont take them
N.sub <- collect(N.long.filt.query, n = Inf) 
N.long.filt.query <- filter(motevo.tbl, gene %in% jgenes.flat)  # peaks are not indexed, so dont take them
N.sub.flat <- collect(N.long.filt.query, n = Inf)

print(head(N.sub))
print(head(N.sub.flat))

print(Sys.time() - start)

stop("STOP")
