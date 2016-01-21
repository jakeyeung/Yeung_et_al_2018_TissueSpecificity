# 2015-12-10
# Jake Yeung
# same as dhs_peak_gene_body.R but across all motifs.

library(hash)
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("/home/yeung/projects/tissue-specificity")
source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
# Load objs -------------------------------------------------------------

load("Robjs/N.long.liver_genes.all_motifs.Robj", verbose=T)
N.long$model <- "Liver"
N.long.all <- N.long
load("Robjs/N.long.flat_genes.all_motifs.Robj", verbose=T)
N.long$model <- "Flat"
N.long.all <- rbind(N.long.all, N.long)
N.long.all <- subset(N.long.all, dist < 1000)

save(N.long.all, file = "Robjs/N.long.flatliver_genes.all_motifs.dist1000.Robj")
