# 2016-03-16
# After running get_tissue_spec_peaks.liver.L.original or KL, plot the UCSC browser of hits in order to understand where 
# the hits are coming from

setwd("/home/yeung/projects/tissue-specificity")
source("scripts/functions/PlotUCSC.R")

# Load --------------------------------------------------------------------

# inf <- "/home/yeung/projects/tissue-specificity/bedfiles/lda_analysis/fg_livrhyth_peaks.bed"
inf <- "/home/yeung/projects/tissue-specificity/bedfiles/lda_analysis/fg_livkidrhyth_peaks.tmp.bed"

peaks.bed <- read.table(inf, col.names = c("chromo", "start", "end", ".", "id"))


# To UCSC -----------------------------------------------------------------

outdir <- "/home/yeung/projects/tissue-specificity/plots/ucsc_motif_screenshots/penalized_lda_peaks"
outname <- "liverkidney_peaks.pdf"
dir.create(outdir)
bedToUCSC(peaks.bed, outpdf = file.path(outdir, outname), leftwindow = 100, rightwindow = 100, theURL = "http://www.genome.ucsc.edu/cgi-bin/hgTracks?hgsid=482393371_BpbPTMseoaFaEGXGPGJTtjOz60ez&hgt.psOutput=on")
