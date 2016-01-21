# 2016-01-07
# Download UCSC browser of liver peaks of Liver-rhythmic and Flat-genes 

setwd("/home/yeung/projects/tissue-specificity")
library(hash)

source("scripts/functions/PlotUCSC.R")

# Functions ---------------------------------------------------------------


CoordToBed <- function(coord){
  chromo <- strsplit(coord, ":")[[1]][[1]]
  startend <- strsplit(coord, ":")[[1]][[2]]
  start <- strsplit(startend, "-")[[1]][[1]]
  end <- strsplit(startend, "-")[[1]][[2]]
  return(data.frame(chromo = chromo, start = as.numeric(start), end = as.numeric(end)))
  # return(matrix(c(chromo, as.numeric(start), as.numeric(end)), ncol = 3, nrow = 1))
}


# Get peaks ---------------------------------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/fits.relamp.Robj", verbose=T)
load("Robjs/S.collapse.liver.dist1000.Robj", verbose=T)
load("Robjs/S.collapse.liverflat.liver_peaks_only.Robj", verbose=T)  # liver peaks for flat genes

# Get liver peaks in liver-rhyth genes
liver.rhyth.livpeaks <- subset(S.collapse, peak.type == "Liver")

# order genes by most rhythmic
key <- as.character(fits.best$gene)
ampval <- as.numeric(fits.best$amp.avg)
ampdic <- hash(key, ampval)

# rename Katnal1;Katnal1 to Katnal because it gives NULL otherwise
liver.rhyth.livpeaks$gene[liver.rhyth.livpeaks$gene == "Katnal1;Katnal1"] <- "Katnal1"
liver.rhyth.livpeaks$amp <- sapply(liver.rhyth.livpeaks$gene, function(g) ampdic[[as.character(g)]], USE.NAMES = FALSE, simplify = TRUE)
liver.rhyth.livpeaks <- liver.rhyth.livpeaks[order(liver.rhyth.livpeaks$amp, decreasing = TRUE), ]

# Get liver peaks in flat-rhyth genes
flat.livpeaks <- subset(S.collapse.livflat, model == "Flat")

# order genes by most expressed in liver
fits.sub <- subset(fits.relamp, tissue == "Liver")
livkey <- as.character(fits.sub$gene)
livval <- as.numeric(fits.sub$int.rnaseq)
livdic <- hash(livkey, livval)

flat.livpeaks$livexprs <- sapply(flat.livpeaks$gene, function(g) {
  exprs <- livdic[[as.character(g)]]
  if (is.null(exprs)){
    return(0)
  } else {
    return(exprs)
  }
})
flat.livpeaks <- flat.livpeaks[order(flat.livpeaks$livexprs, decreasing = TRUE), ]

# Do top n genes ---------------------------------------------------------

top.n <- 50
genes.sub <- head(unique(liver.rhyth.livpeaks$gene), n = top.n)
genes.flat.sub <- head(unique(flat.livpeaks$gene), n = top.n)

# for liver rhyth genes
liver.rhyth.livpeaks.bed <- lapply(as.character(subset(liver.rhyth.livpeaks, gene %in% genes.sub)$peak), function(peak){
  CoordToBed(as.character(peak))
})
liver.rhyth.livpeaks.bed <- do.call(rbind, liver.rhyth.livpeaks.bed)
head(liver.rhyth.livpeaks.bed)
str(liver.rhyth.livpeaks.bed)

# for flat genes
flat.livpeaks.bed <- lapply(as.character(subset(flat.livpeaks, gene %in% genes.flat.sub)$peak), function(peak){
  CoordToBed(as.character(peak))
})
flat.livpeaks.bed <- do.call(rbind, flat.livpeaks.bed)


# Plot UCSC ---------------------------------------------------------------

bedToUCSC(toPlot = liver.rhyth.livpeaks.bed, 
          outpdf = "/home/yeung/projects/tissue-specificity/plots/ucsc_motif_screenshots/motif_views_liver_peaks_liver_rhyth.pdf", 
          leftwindow = 100, rightwindow = 100, 
          theURL = "http://www.genome.ucsc.edu/cgi-bin/hgTracks?hgsid=467519929_OQqiHWdViYPcoVoWUN8mIPEy41H5&hgt.psOutput=on",
          chromo.name="chromo", start.name="start", end.name="end")

bedToUCSC(toPlot = flat.livpeaks.bed, 
          outpdf = "/home/yeung/projects/tissue-specificity/plots/ucsc_motif_screenshots/motif_views_liver_peaks_flat.pdf", 
          leftwindow = 100, rightwindow = 100, 
          theURL = "http://www.genome.ucsc.edu/cgi-bin/hgTracks?hgsid=467519929_OQqiHWdViYPcoVoWUN8mIPEy41H5&hgt.psOutput=on",
          chromo.name="chromo", start.name="start", end.name="end")
