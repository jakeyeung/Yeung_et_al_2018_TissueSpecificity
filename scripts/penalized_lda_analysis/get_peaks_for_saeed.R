# 2016-02-19
# Saeed wants peaks to be compared between two conditions
# do it in bedformat

library(reshape2)

# Functions ---------------------------------------------------------------

PeaksToBed <- function(peaks, name){
  # Create bed file given peaks and peak name
  colsplit <- lapply(peaks, function(p) strsplit(p, ":")[[1]])
  chromo <- sapply(colsplit, function(jsplit) jsplit[[1]])
  start <- sapply(colsplit, function(jsplit) strsplit(jsplit[[2]], "-")[[1]][[1]])
  end <- sapply(colsplit, function(jsplit) strsplit(jsplit[[2]], "-")[[1]][[2]])
  dat <- data.frame(chromo = chromo, 
                    start = start,
                    end = end,
                    id = name)
  return(dat)
}


# Load --------------------------------------------------------------------

# do for peaksets

indir <- "/home/shared/tissue_specificity/peak_sets"
inf <- "S.liver.liverspec.peaks.minamp02.Robj"
inf <- "S.flat.liverspec.peaks.top70percent.Robj"
bname <- strsplit(inf, ".Robj")[[1]][[1]]
outf <- paste0(bname, ".txt")
load(file.path(indir, inf), verbose=T)

bed.liverspec <- PeaksToBed(as.character(S.sub.collapse.filt$peak), as.character(S.sub.collapse.filt$gene))

write.table(bed.liverspec, file = file.path(indir, outf), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# do for sitecounts set

inf.sc <- "/home/shared/tissue_specificity/sitecounts_sets/N.long.sum.bytiss.Robj"
outf.sc <- "/home/shared/tissue_specificity/sitecounts_sets/N.long.sum.bytiss.txt"
load(inf.sc, v=T)  # N.long.sum.bytiss

N.mat <- dcast(N.long.sum.bytiss, formula = gene ~ tissue + motif, value.var = "sitecount.max", fill = 0)

write.table(N.mat, file = outf.sc, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
