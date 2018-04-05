# 2016-02-26
# Jake Yeung
# We found RXRG_dimer.p3, HNF1A, FOXA2.p3, ONECUT1,2.p2 in our hits when we summed motifs across gene body
# What kind of peaks are in RXRG, HNF1A, FOXA2, ONECUT1? 

library(hash)
library(dplyr)
library(ggplot2)

setwd("/home/yeung/projects/tissue-specificity")

# Functions ---------------------------------------------------------------

source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/DataHandlingFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotUCSC.R")

source("scripts/functions/LdaFunctions.R")
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
source("scripts/functions/PlotGeneAcrossTissues.R")

load("Robjs/N.long.sum.bytiss.all_genes.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)

library(reshape2)

# Load --------------------------------------------------------------------

load("Robjs/S.long.Robj", verbose=T)
load("Robjs/N.long.all_genes.all_signif_motifs.Robj", v=T)
N.long.liver_dhs <- N.long.filt; rm(N.long.filt)
load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)


# Run LDA -----------------------------------------------------------------

# Load --------------------------------------------------------------------

minamp <- 0.25
sitecount.name <- "sitecount.max"
bg.genes.type <- "Flat"  # "Flat", "Rhythmic", "same" 


# Kidney Liver models
jgrep <- "^Kidney;Liver$|^Kidney,Liver$"
jmodels <- unique(fits.best[grepl(jgrep, fits.best$model), ]$model)

# define foreground and background tissues
all.tiss <- c('Cere', 'Heart', 'Kidney', 'Liver', 'Lung', 'Mus')
fg.tiss <- c("Liver", "Kidney")
bg.tiss <- all.tiss[which(! all.tiss %in% fg.tiss & all.tiss != "Kidney")]

out <- RunPenalizedLDA(N.long.sum.bytiss, fits.best, jmodels, fg.tiss, fg.tiss, bg.genes.type, sitecount.name, minamp = 0.25, n.bg = "all")

cnames <- colnames(out$x)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]

# Sort peaks by presence of jmotif ----------------------------------------

# jmotifs <- c('RXRG_dimer.p3', 'HNF1A.p2', 'FOXA2.p3', 'ONECUT1,2.p2', 'HNF4A_NR2F1,2.p2', 
#              'ESRRA.p2', 'CUX2.p2', 'NR4A2.p2', 'FOX{F1,F2,J1}.p2', 'NR6A1.p2', 'FOXQ1.p2', 
#              'HLF.p2', 'NR5A1,2.p2', 'CDC5L.p2', 'NR1H4.p2', 'RORA.p2', 'SRY.p2', 'FOXO1,3,4.p2', 'ATF2.p2', 'bHLH_family.p2')
# jmotifs <- c('HNF1A.p2','HNF4A_NR2F1,2.p2','ADNP_IRX_SIX_ZHX.p2','EVI1.p2','ESRRA.p2','NR6A1.p2','ATF2.p2','FOX{I1,J2}.p2',
#              'NR5A1,2.p2','FOXD3.p2','FOX{F1,F2,J1}.p2','RUNX1..3.p2','FOXA2.p3','ONECUT1,2.p2','NR1H4.p2','CDX1,2,4.p2','LHX3,4.p2','CUX2.p2','FOXQ1.p2','LEF1_TCF7_TCF7L1,2.p2')
jmotifs <- cnames[1:10]
jtop.n <- 10
# jmotif <- "RXRG_dimer.p3"
# jmotif <- "ONECUT1,2.p2"
# jmotif <- "HNF1A.p2"
# jmotif <- "FOXA2.p3"

start <- Sys.time()
for (jmotif in jmotifs){
  outf <- paste0("/home/yeung/projects/tissue-specificity/plots/ucsc_motif_screenshots/motif_views.kidneyliv.", jmotif, ".pdf")
  print(jmotif)
  if (file.exists(outf)){
    print(paste("Skipping motif", jmotif))
    next 
  } 

  N.long.liver_dhs.liv.motif <- subset(N.long.liver_dhs.liv, motif == jmotif)
  N.long.liver_dhs.liv.motif <- N.long.liver_dhs.liv.motif[order(N.long.liver_dhs.liv.motif$sitecount, decreasing = TRUE), ]
  
  print(nrow(N.long.liver_dhs.liv.motif))
  
  top.n <- min(jtop.n, nrow(N.long.liver_dhs.liv.motif))
  
  N.long.liver_dhs.liv.motif <- N.long.liver_dhs.liv.motif[1:top.n, ]
  
  # Plot peaks  -------------------------------------------------------------
  
  if (nrow(N.long.liver_dhs.liv.motif) > 0){
    # for liver rhyth genes
    peaks.bed <- lapply(as.character(N.long.liver_dhs.liv.motif$peak), function(peak){
      CoordToBed(as.character(peak))
    })
    peaks.bed <- do.call(rbind, peaks.bed)
    
    print(head(peaks.bed))
    bedToUCSC(toPlot = peaks.bed, 
              # outpdf = "/home/yeung/projects/tissue-specificity/plots/ucsc_motif_screenshots/motif_views_rxrg_dimer_peaks.pdf", 
              outpdf = outf, 
              leftwindow = 250, rightwindow = 250, 
              theURL = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=477175389_DSFdGfM7dyhhHORzyaaej8ICXXQQ&hgt.psOutput=on", 
              chromo.name="chromo", start.name="start", end.name="end")
  } else {
    print(paste("Skipping", jmotif, "empty data frame"))
  }
  
}
print(Sys.time() - start)
