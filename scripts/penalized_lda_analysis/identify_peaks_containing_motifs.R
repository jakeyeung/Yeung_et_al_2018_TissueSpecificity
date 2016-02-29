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

# Load --------------------------------------------------------------------

load("Robjs/S.long.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)
load("Robjs/N.long.liver_genes.all_motifs.100000.Robj", verbose=T)  # huge file 
load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)

# Get liver genes ---------------------------------------------------------

ampmin <- 0.25
liv.genes <- unique(as.character(subset(fits.best, model == "Liver" & amp.avg > ampmin)$gene))

# shift by 1 log2 unit the fit looks better
log.shift <- 2.5
S.tissuecutoff$cutoff.adj <- 2^(log2(S.tissuecutoff$cutoff) + log.shift)

cutoff.adj <- S.tissuecutoff$cutoff.adj
cutoff.lower <- S.tissuecutoff$cutoff

# show upper and lower limit
pseudo <- 1e-3
ggplot(S.long[sample(x = 1:nrow(S.long), size = 0.01 * nrow(S.long)), ], aes(x = log2(signal + pseudo))) + geom_density() + 
  facet_wrap(~tissue) +
  geom_vline(aes(xintercept = log2(cutoff.adj)), data = S.tissuecutoff, colour = "blue") + 
  geom_vline(aes(xintercept = log2(cutoff)), data = S.tissuecutoff, colour = "red")

cutoffs.tiss.upper <- hash(as.character(S.tissuecutoff$tissue), as.numeric(S.tissuecutoff$cutoff.adj))

jstart <- Sys.time()  # 4 min 
S.sub.collapse <- subset(S.long, gene %in% liv.genes)

# apply is faster than mapply
S.sub.collapse$is.upper <- apply(S.sub.collapse, 1, function(row) IsSignalUpper(as.character(row[5]), as.numeric(row[6]), cutoffs.tiss.upper))

S.sub.collapse.peaks <- subset(S.sub.collapse, is.upper == TRUE)

(nrow(S.sub.collapse))
(nrow(S.sub.collapse.peaks))


# Filter N.long for peaks  ------------------------------------------------

# need to do it tissue by tissue
tissues <- as.character(unique(S.sub.collapse.peaks$tissue))

peaks.bytissue <- lapply(tissues, function(tiss){
  unique(as.character(subset(S.sub.collapse.peaks, tissue == tiss)$peak))
})
names(peaks.bytissue) <- tissues


# Find peaks in Liver  ----------------------------------------------------

jtiss <- "Liver"

N.long.liver_dhs.liv <- subset(N.long.liver_dhs, peak %in% peaks.bytissue[[jtiss]])


# Sort peaks by presence of jmotif ----------------------------------------

jmotifs <- c('RXRG_dimer.p3', 'HNF1A.p2', 'FOXA2.p3', 'ONECUT1,2.p2', 'HNF4A_NR2F1,2.p2', 
             'ESRRA.p2', 'CUX2.p2', 'NR4A2.p2', 'FOX{F1,F2,J1}.p2', 'NR6A1.p2', 'FOXQ1.p2', 
             'HLF.p2', 'NR5A1,2.p2', 'CDC5L.p2', 'NR1H4.p2', 'RORA.p2', 'SRY.p2', 'FOXO1,3,4.p2', 'ATF2.p2', 'bHLH_family.p2')
# jmotif <- "RXRG_dimer.p3"
# jmotif <- "ONECUT1,2.p2"
# jmotif <- "HNF1A.p2"
# jmotif <- "FOXA2.p3"

start <- Sys.time()
for (jmotif in jmotifs){
  outf <- paste0("/home/yeung/projects/tissue-specificity/plots/ucsc_motif_screenshots/motif_views.", jmotif, ".pdf")
  print(jmotif)
  if (file.exists(outf)){
    print(paste("Skipping motif", jmotif))
    next 
  } 

  N.long.liver_dhs.liv.motif <- subset(N.long.liver_dhs.liv, motif == jmotif)
  N.long.liver_dhs.liv.motif <- N.long.liver_dhs.liv.motif[order(N.long.liver_dhs.liv.motif$sitecount, decreasing = TRUE), ]
  
  print(nrow(N.long.liver_dhs.liv.motif))
  
  top.n <- min(100, nrow(N.long.liver_dhs.liv.motif))
  
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
