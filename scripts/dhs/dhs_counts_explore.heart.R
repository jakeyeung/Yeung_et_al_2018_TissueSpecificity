# Look at DHS counts data
# Date: 2015-03-17
rm(list=ls())

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


dhs.path <- "/home/yeung/projects/tissue-specificity/data/beds/ENCODE:Seq:Heart/combined.bed"
dhs.dat <- read.table(dhs.path, header = TRUE, sep = "\t")  # gigantic file ~ 1.2 gigs in memory
# save(dhs.dat, file = "~/projects/tissue-specificity/docs/2015-03-19-felix/dhs_dat.Robj")
jtissue <- "Heart"

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


# Remove replicate 1, keep replicate 2 and 3 ------------------------------

outlier <- "UwStam_mHeart.DS18138.FC6323D.1_001"

good.samples <- colnames(dhs.reps.sub)[which(colnames(dhs.reps.sub) != outlier)]

head(dhs.reps[, good.samples])


# Find cutoff -------------------------------------------------------------

source("scripts/functions/MixtureModelFunctions.R")

counts <- unlist(c(dhs.reps[rows.sub, good.samples]))
counts <- log2(as.numeric(counts[which(counts > 0)]))
plot(density(counts))

cutoff.log2 <- FindCutoff(x = counts, lambdas = c(0.6, 0.4), mus = c(-5, 0))
abline(v = cutoff.log2$maximum)  # should intersect two gaussians

cutoff <- 2^cutoff.log2$maximum

print(paste("Cutoff:", cutoff, "log2 Cutoff:", cutoff.log2$maximum))

# Take mean of good samples, filter for cutoff ----------------------------

dhs.clean <- data.frame(chr = dhs.dat$chr,
                        start = dhs.dat$start,
                        end = dhs.dat$end,
                        filtered_norm_counts = apply(dhs.reps[, good.samples], 1, mean))

dhs.clean.filtered <- dhs.clean[which(dhs.clean$filtered_norm_counts > cutoff), ]

plot(hist(log2(dhs.clean.filtered$filtered_norm_counts), 100))

# Write to output ---------------------------------------------------------

write.table(dhs.clean.filtered, file = "data/beds/filtered_beds/encode_peaks.heart.bed", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
