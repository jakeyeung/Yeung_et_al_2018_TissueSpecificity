# 2015-12-09
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


inf <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_stringent/merged/merged.all.bed"
cnames <- c("chromo", "start", "end", "signal", "tissue")
S.stringent <- read.table(inf, header = FALSE, sep = "\t", col.names = cnames)

# make chromo,start,end a single column
# mapply is slow
S.stringent$peak <- mapply(MakeCoord, S.stringent$chromo, S.stringent$start, S.stringent$end)

print(head(S.stringent))
save(S.stringent, file = "Robjs/S.stringent.Robj")

print(Sys.time() - start)