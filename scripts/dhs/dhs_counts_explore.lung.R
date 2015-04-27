# Look at DHS counts data
# Date: 2015-03-17
# rm(list=ls())

library(ggplot2)

# Functions ---------------------------------------------------------------

GetReplicates <- function(cnames, get.index=TRUE){
  # from colnames of data, retrieve replicate names
  # expects that colnames are not called "chr", "start", "end", "zscore" or "total"
  # get.index -> return the index from cnames otherwise get names
  colnames.remove <- c("chr", "start", "end", "zscore", "total")
  if (get.index){
    return(which(!cnames %in% colnames.remove))
  }
  else {
    return(cnames[which(!cnames %in% colnames.remove)])
  }
}

fit <- function(x, a, b){
  return(exp(a*x + b))
}

# Main --------------------------------------------------------------------


# Get data ----------------------------------------------------------------


dhs.path <- "/home/yeung/projects/tissue-specificity/data/beds/ENCODE:Seq:Lung/combined.bed"
dhs.dat <- read.table(dhs.path, header = TRUE, sep = "\t")  # gigantic file ~ 1.2 gigs in memory
# save(dhs.dat, file = "~/projects/tissue-specificity/docs/2015-03-19-felix/dhs_dat.Robj")
jtissue <- "Lung"

# subset constants --------------------------------------------------------
set.seed(0)
n.sub <- 0.01 * nrow(dhs.dat)
rows.sub <- sample(seq(1:nrow(dhs.dat)), n.sub)
i.reps <- GetReplicates(colnames(dhs.dat))
n.samps <- length(i.reps)

# Normalize by total counts -----------------------------------------------

dhs.reps.sub <- dhs.dat[rows.sub, i.reps]
sum.samps <- apply(dhs.dat[, i.reps], 2, sum)
dhs.reps <- dhs.dat[, i.reps] / sum.samps * 10^6

barplot(sum.samps, las = 2)

pairs(log2(dhs.reps.sub))

for (i in 1:ncol(dhs.reps.sub)){
  plot(density(log2(unlist(dhs.reps.sub[, i]))), main = paste(jtissue, "rep", i))
}

