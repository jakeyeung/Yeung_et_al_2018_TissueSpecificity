# 2015-09-06
# Jake Yeung
# dhs_filter_by_gamma.R


# Functions ---------------------------------------------------------------

source("scripts/functions/MixtureModelFunctions.R")

GetBedPath <- function(tissue, beddir="/home/yeung/projects/tissue-specificity/data/beds/merge/ucsc_names"){
  bedfname <- paste(tissue, "dhs", "merged", "bed", sep = ".")
  bedpath <- file.path(beddir, bedfname)
  return(bedpath)
}

# Load data ---------------------------------------------------------------


tissues <- c("SkeletalMuscle", "Liver", "Cerebellum", "Lung", "Kidney", "Heart")

for (tissue in tissues){
  print(tissue)
  beddir <- "/home/yeung/projects/tissue-specificity/data/beds/merge/ucsc_names"
  bedpath <- GetBedPath(tissue, beddir)
  
  gammadir <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_smoother_gamma_cutoff_high_accuracy_ALL_gamma_inits"
  fname <- "gammamdl.Robj"
  gammapath <- file.path(gammadir, tissue, fname)
  
  outfname <- paste(tissue, "gamma_filt", "bed", sep = ".")
  outbed <- file.path(gammadir, tissue, outfname)
  
  dat <- read.table(bedpath, header = FALSE)
  
  counts <- dat$V4[which(dat$V4 > 1)]
  counts <- sample(counts, 0.01 * length(counts))  
  
  load(gammapath, verbose = T)
  
  # Get cutoff --------------------------------------------------------------
  
  cutoff.log2 <- PlotGammaMixmdl(mixmdl, log2(counts), savedir = file.path(gammadir, tissue))
  
  cutoff <- 2 ^ cutoff.log2
  
  dat.filt <- dat[which(dat$V4 >= cutoff), ]
  write.table(dat.filt, file = outbed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

