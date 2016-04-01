# 2016-03-18
# Convert signal to Z-score (more homogeneous)

rm(list=ls())
library(dplyr)
library(ggplot2)


# Functions ---------------------------------------------------------------

SignalToBed <- function(jsub, outdir){
  jtiss <- as.character(jsub$tissue[[1]])
  # change tissues to match RNA-Seq names
  if (jtiss == "Mus") jtiss <- "SkeletalMuscle"
  if (jtiss == "Cere") jtiss <- "Cerebellum"
  dir.create(outdir, showWarnings = FALSE)
  outfile <- file.path(outdir, paste0(jtiss, ".dhs.noMT.bed"))  # noMT because it fits with downstream scripts
  write.table(subset(jsub, select = c(chromo, start, end, signal.std)), file = outfile, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  return(data.frame())
}

# Load --------------------------------------------------------------------

load("Robjs/S.long.Robj", v=T)


# Standardize -------------------------------------------------------------

pseudo <- 0.01
S.long <- S.long %>%
  group_by(tissue) %>%
  mutate(signal.std = scale(log2(signal + pseudo), center = TRUE, scale = TRUE))


# Plot to check ----------------------------------------------------------

ggplot(S.long[sample(seq(nrow(S.long)), size = 0.01 * nrow(S.long)), ], aes(x = signal.std)) + 
  geom_density() + facet_wrap(~tissue) + theme_bw(24) + geom_vline(aes(xintercept=0)) + geom_hline(aes(yintercept=0))


# Add cutoff --------------------------------------------------------------

# be stringent

cutoff <- 1.8
ggplot(S.long[sample(seq(nrow(S.long)), size = 0.01 * nrow(S.long)), ], aes(x = signal.std)) + 
  geom_density() + facet_wrap(~tissue) + theme_bw(24) + geom_vline(aes(xintercept=cutoff)) + geom_hline(aes(yintercept=0))



# Print peaks with cutoffs ------------------------------------------------

outdir <- "/home/yeung/projects/tissue-specificity/data/beds/merge/cutoffs_stringent"

out <- S.long %>%
  filter(signal.std > cutoff) %>%
  group_by(tissue) %>%
  do(SignalToBed(., outdir))

