# Differential sitecount analysis: different within tissue-specific modules
# 2015-11-20

setwd("~/projects/tissue-specificity/")

library(hash)
library(dplyr)
library(ggplot2)
library(reshape2)

source("scripts/functions/LoadSitecounts.R")
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ReadListToVector.R")
source("scripts/functions/DifferentialSitecountsFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")

# Functions ---------------------------------------------------------------



# Load --------------------------------------------------------------------


# N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_matrix"
# suffix <- "sitecounts.merged.matrix"
# N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_mean_matrix"
# suffix <- "sitecounts.dist_filtered.mean.matrix"
# N <- LoadSitecountsEncodeAll(maindir = N.dir, suffix = suffix, with.ensemblid = FALSE, rename.tissues = TRUE)  # merged by gene
# N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_mean_matrix_bugfixed_redo"
N.dir <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/encode_dist_filtered_sum_matrix_bugfixed"
suffix <- "dist_filt.bugfixed.sitecount.matrix"
N <- LoadSitecountsEncodeAll(maindir = N.dir, tissues = c("Liver", "Kidney", "Cere", "Lung", "Heart", "Mus"),
                             suffix = suffix, with.ensemblid = FALSE, rename.tissues = FALSE)  # merged by gene


load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj", verbose=T)

genes.liver <- as.character(subset(fits.best, model == "Liver")$gene)
genes.flat <- as.character(subset(fits.best, model == "")$gene)


# Do 2 by 2 tables ---------------------------------------------------------

N.sub <- subset(N, gene %in% genes.liver)

motifs <- as.character(unique(N.sub$motif))

N.sub$sitecount <- N.sub$motevo.value
N.sub$model <- sapply(as.character(N.sub$tissue), function(x){
  if (x == "Liver"){
    return("Liver")
  } else {
    return("NotLiver")
  }
})

test.out <- N.sub %>%
  group_by(motif) %>%
  do(FisherTestSitecounts(., cutoff = 0.6))
test.out[order(test.out$p.value), ]


start <- Sys.time()
cutoffs <- seq(from = 0.5, to = 2, by = 0.5)
test.all <- data.frame()
for (cutoff in cutoffs){
  print(cutoff)
  test.out <- N.sub %>%
    group_by(motif) %>%
    do(FisherTestSitecounts(., cutoff))
  test.out$cutoff <- cutoff
  test.all <- rbind(test.all, test.out)
}
print(Sys.time() - start)

test.sum <- test.all %>%
  group_by(motif) %>%
  summarise(odds.ratio = mean(odds.ratio), p.value.minuslog = mean(-log10(p.value)))

ggplot(test.sum, aes(x = p.value.minuslog, y = odds.ratio, label = motif)) + geom_point() + geom_text()

FisherTestSitecounts(subset(N.sub, motif == "CUX2.p2"), cutoff=1, show.table = TRUE)
