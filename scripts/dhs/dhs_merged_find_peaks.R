# dhs_merged_find_peaks.R
# 2015-05-28
# Jake Yeung

# Get parameters ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
fname <- file.path(args[1]) 
# fname <- "data/beds/merge/ucsc_names/Liver.dhs.merged.bed"
outdir <- file.path(args[2])  # no extension, add it in later
# outdir <- "data/beds/merge/Liver.cutoff"
outbed <- file.path(args[3])  # full path yo
filter_counts <- as.numeric(args[4])
init_mu1 <- as.numeric(args[5])
init_mu2 <- as.numeric(args[6])

fname.split <- strsplit(basename(fname), "\\.")[[1]]
fname.base <- paste(fname.split[1:(length(fname.split) - 1)], collapse = ".")

print(fname.base)
print(paste0(outdir, '/', fname.base, ".pdf"))

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

counts <- dat$V4[which(dat$V4 > 1)]
counts <- sample(counts, 0.01 * length(counts))  # make curve a little smoother


# outdir <- "data/beds/merge/Liver.cutoff.pdf"
print("Begin: fit mixture model")
start <- Sys.time()
cutoff.log2 <- FindCutoff(x = log2(counts), lambdas = c(0.7, 0.3), mus = c(init_mu1, init_mu2), outdir=outdir)
# fit gammae
gscale <- 50
gmeans <- c(3, 12)
glambdas <- c(0.8, 0.2)
mixmdl <- gammamixEM(log2(counts), lambda = glambdas, alpha = c(gscale, gscale), beta = c(gmeans[1] / gscale, gmeans[2] / gscale))
save(mixmdl, file = file.path(outdir, "gammamdl.Robj"))
# find cutoff for gamma
cutoff.gamma.log2 <- PlotGammaMixmdl(mixmdl, log2(counts), savedir = outdir)
cutoff.gamma <- 2 ^ cutoff.gamma.log2
print("End: fit mixture model")
print(Sys.time() - start)

pdf(paste0(outdir, "/", fname.base, ".pdf"))
plot(density(log2(counts)))
abline(v = cutoff.log2$maximum)
plot(mixmdl, which = 2)
dev.off()

cutoff <- 2^cutoff.log2$maximum
print(paste("Cutoff:", cutoff, "log2 Cutoff:", cutoff.log2$maximum))


# Filter by cutoff --------------------------------------------------------

dat.filt <- dat[which(dat$V4 >= cutoff), ]
write.table(dat.filt, file = outbed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Cutoff by making all data that do not pass cutoffs to zero. This makes merging easier.
dat.filt <- dat
dat.filt$V4[which(dat.filt$V4 < cutoff.gamma)] <- 0
write.table(dat.filt, file = file.path(outdir, paste0(fname.base, ".gamma_filtered.bed")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Normalize by total counts
norm.factor <- sum(dat$V4) / 10^9
dat.filt$V4 <- dat.filt$V4 / norm.factor
write.table(dat.filt, file = file.path(outdir, paste0(fname.base, ".gamma_filtered.normalized.bed")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Plot filtered distributions ---------------------------------------------

pdf(file.path(outdir, paste0(fname.base, ".filtered_density.pdf")))
plot(density(log2(dat.filt$V4)))
dev.off()
