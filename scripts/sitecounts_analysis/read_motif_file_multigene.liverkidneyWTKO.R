# 2016-07-08
# Jake Yeung
# read_filter_motif_files.R
# read each motif file and filter it for gene list, then concatenate into a larger file and save to output

start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")


# Functions ---------------------------------------------------------------

MakeCoord <- function(chromo, start, end){
  startend = paste0(start, "-", end)
  coord = paste0(chromo, ":", startend)
  return(coord)
}

# Main --------------------------------------------------------------------



inf <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed_multigene/all_sites.closest.filter.0.5.dist.25000.long.bed"
cnames <- c("chromo", "start", "end", "motif", "sitecount", "chromo.peak", "start.peak", "end.peak", "gene", "dist", "")
N.multigene <- read.table(inf, header = FALSE, sep = "\t", col.names = cnames)
N.multigene$X <- NULL

# make chromo,start,end a single column
# mapply is slow
N.multigene$coord.motif <- mapply(MakeCoord, N.multigene$chromo, N.multigene$start, N.multigene$end)
N.multigene$coord.peak <- mapply(MakeCoord, N.multigene$chromo.peak, N.multigene$start.peak, N.multigene$end.peak)

save(N.multigene, file = "Robjs/N.multigene.distfilt.25000.Robj")

print(Sys.time() - start)
