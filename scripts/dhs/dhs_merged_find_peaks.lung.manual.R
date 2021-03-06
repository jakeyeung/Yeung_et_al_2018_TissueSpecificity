# dhs_merged_find_peaks.R
# 2015-05-28
# Jake Yeung
rm(list=ls())
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

# fname <- "data/beds/merge/ucsc_names/Heart.dhs.merged.bed"
fname <- "data/beds/merge/ucsc_names/Lung.dhs.merged.bed"
outbed <- "data/beds/merge/cutoffs_manual/Lung.dhs.merged.filtered.manual.bed"
dat <- read.table(fname, header = FALSE)

# Find peaks from merged bam files ----------------------------------------

counts <- dat$V4[which(dat$V4 > 0)]
counts <- sample(counts, length(counts))
plot(density(log2(counts)))

cutoff.log2 <- 5.9
abline(v = cutoff.log2)
cutoff <- 2^cutoff.log2
print(paste("Cutoff:", cutoff))

# # outdir <- "data/beds/merge/Liver.cutoff.pdf"
# print("Begin: fit mixture model")
# start <- Sys.time()
# # mixmdl = normalmixEM(log2(counts), lambda = c(0.3, 0.3, 0.3), mu = c(3, 5, 7), k = 3)
# # cutoff <- optimize(f = ShannonEntropyMixMdl, interval = range(x), mixmdl = mixmdl, tol = 0.0001, maximum = TRUE)
# cutoff.log2 <- FindCutoff(x = log2(counts), lambdas = c(0.5, 0.5), mus = c(4, 6))
# # cutoff.log2 <- FindCutoff(x = log2(counts), lambdas = c(1, 1, 0.1), mus = c(3, 7, 5), k = 3)
# print("End: fit mixture model")
# print(Sys.time() - start)
# cutoff <- 2^cutoff.log2$maximum
# print(paste("Cutoff:", cutoff, "log2 Cutoff:", cutoff.log2$maximum))


# Filter by cutoff --------------------------------------------------------

dat.filt <- dat[which(dat$V4 >= cutoff), ]
# remove chrMT
dat.filt <- dat.filt[which(dat.filt$V1 != "chrMT"), ]
write.table(dat.filt, file = outbed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

