# Jake Yeung
# 2015-09-03 
# Test out Gamma fitting functions

source("scripts/functions/MixtureModelFunctions.R")

bedpath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/ucsc_names/Cerebellum.dhs.merged.bed"
gammapath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma/Cerebellum/gammamdl.Robj"
gammapath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy/Cerebellum/gammamdl.Robj"
gammapath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma_inits/Cerebellum/gammamdl.Robj"

bedpath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/ucsc_names/Heart.dhs.merged.bed"
gammapath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma2/Heart/gammamdl.Robj"
gammapath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma_inits/Heart/gammamdl.Robj"

bedpath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/ucsc_names/SkeletalMuscle.dhs.merged.bed"
gammapath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma_inits/SkeletalMuscle/gammamdl.Robj"

bedpath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/ucsc_names/Liver.dhs.merged.bed"
gammapath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma_inits/Liver/gammamdl.Robj"
gammapath <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma_inits_scale50/Liver/gammamdl.Robj"

dat <- read.table(bedpath, header = FALSE)

# Find peaks from merged bam files -------- -------------------------------

counts <- dat$V4[which(dat$V4 > 1)]
counts <- sample(counts, 0.01 * length(counts))


# Load object -------------------------------------------------------------

load(gammapath, verbose = TRUE)

# PlotGammaDist(mixmdl, log2(counts))

# scale <- 50
# mixmdl <- gammamixEM(log2(counts), lambda = c(0.8, 0.2), alpha = c(scale, scale), beta = c(3 / scale, 12 / scale))
cutoff <- PlotGammaMixmdl(mixmdl, log2(counts))

# cutoff <- optimize(f = ShannonEntropyMixMdl, c(log2(counts), mixmdl), interval = range(log2(counts)), tol = 0.0001, maximum = TRUE)

# abline(v = cutoff$maximum)
