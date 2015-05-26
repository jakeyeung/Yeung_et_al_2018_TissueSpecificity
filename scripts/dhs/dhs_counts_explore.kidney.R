# Look at DHS counts data
# Date: 2015-03-17
rm(list=ls())

library(ggplot2)

# Functions ---------------------------------------------------------------

source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/DhsFunctions.R")

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


# dhs.path <- "~/projects/tissue-specificity/data/beds/ENCODE:Seq:Liver/combined.nogc.bedfixed"
dhs.path <- "/home/yeung/projects/tissue-specificity/data/beds/ENCODE:Seq:Kidney/combined.bed"
dhs.dat <- read.table(dhs.path, header = TRUE, sep = "\t")  # gigantic file ~ 1.2 gigs in memory
# save(dhs.dat, file = "~/projects/tissue-specificity/docs/2015-03-19-felix/dhs_dat.Robj")
jtissue <- "Kidney"

# subset constants --------------------------------------------------------
set.seed(0)
n.sub <- 0.01 * nrow(dhs.dat)
rows.sub <- sample(seq(1:nrow(dhs.dat)), n.sub)
i.reps <- GetReplicates(colnames(dhs.dat))
n.samps <- length(i.reps)

# Normalize by total counts -----------------------------------------------

sum.samps <- colSums(dhs.dat[, i.reps])
dhs.reps <- sweep(dhs.dat[, i.reps], 2, sum.samps, "/") * 10^6

barplot(sum.samps, las = 2)

dhs.reps.sub <- dhs.dat[rows.sub, i.reps]  # for plotting

pairs(log2(dhs.reps.sub))

for (i in 1:ncol(dhs.reps.sub)){
  print(colnames(dhs.reps.sub)[i])
  plot(density(log2(unlist(dhs.reps.sub[, i]))), main = paste(jtissue, colnames(dhs.reps.sub)[i]))
}


# Find cutoff -------------------------------------------------------------

good.samples <- colnames(dhs.reps.sub)  # all good
counts <- unlist(c(dhs.reps[rows.sub, good.samples]))
counts <- log2(as.numeric(counts[which(counts > 0)]))
plot(density(counts))

cutoff.log2 <- FindCutoff(x = counts, lambdas = c(0.6, 0.4), mus = c(-5, -1))
cutoff <- 2^cutoff.log2$maximum
print(paste("Cutoff:", cutoff, "log2 Cutoff:", cutoff.log2$maximum))


# Filter output by cutoff (take mean of replicates) -----------------------

dhs.clean.filtered <- FilterReadcounts(dhs.dat, dhs.reps, good.samples, cutoff)

# Write to file -----------------------------------------------------------

write.table(dhs.clean.filtered, file = "data/beds/filtered_beds/encode_peaks.kidney.bed", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
