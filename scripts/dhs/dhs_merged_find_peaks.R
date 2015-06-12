# dhs_merged_find_peaks.R
# 2015-05-28
# Jake Yeung

# Get parameters ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
fname <- file.path(args[1]) 
# fname <- "data/beds/merge/ucsc_names/Liver.dhs.merged.bed"
outdir <- file.path(args[2])  # no extension, add it in later
# outdir <- "data/beds/merge/Liver.cutoff"
outbed <- file.path(args[3])
init_mu1 <- as.numeric(args[4])
init_mu2 <- as.numeric(args[5])

if (is.na(outbed)){
  outbed <- file.path(outdir, "filtered.bed")
  sprintf("No bed specified, defaulting: %s", outbed)
}
if (is.na(init_mu1) | is.na(init_mu2)){
  init_mu1 <- 3
  init_mu2 <- 7
  sprintf("No init mu specified, using %s, %s", init_mu1, init_mu2)
}

print(args)


# Create output directory -------------------------------------------------

dir.create(outdir, showWarnings = FALSE)

# Set directory and functions ---------------------------------------------


setwd("~/projects/tissue-specificity")
library(mixtools)
source("scripts/functions/MixtureModelFunctions.R")

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}


# Load file ---------------------------------------------------------------

# fname <- "data/beds/merge/ucsc_names/Liver.dhs.merged.bed"
dat <- read.table(fname, header = FALSE)

# Find peaks from merged bam files ----------------------------------------

counts <- dat$V4[which(dat$V4 > 0)]
counts <- sample(counts, 0.01 * length(counts))  # make curve a little smoother

plot(density(log2(counts)))

# outdir <- "data/beds/merge/Liver.cutoff.pdf"
print("Begin: fit mixture model")
start <- Sys.time()
cutoff.log2 <- FindCutoff(x = log2(counts), lambdas = c(0.7, 0.3), mus = c(init_mu1, init_mu2), outdir=outdir)
print("End: fit mixture model")
print(Sys.time() - start)


cutoff <- 2^cutoff.log2$maximum
print(paste("Cutoff:", cutoff, "log2 Cutoff:", cutoff.log2$maximum))


# Filter by cutoff --------------------------------------------------------

dat.filt <- dat[which(dat$V4 >= cutoff), ]
write.table(dat.filt, file = outbed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

