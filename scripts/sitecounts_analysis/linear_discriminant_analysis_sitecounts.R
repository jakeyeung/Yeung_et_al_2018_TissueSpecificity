# 2015-12-17
# Find combination of motifs that separate between my groups

library(hash)
library(MASS)
library(reshape2)
source("scripts/functions/FisherTestSitecounts.R")

# Load --------------------------------------------------------------------

load("Robjs/N.long.liver_genes.all_motifs.Robj", verbose=T)
load("Robjs/S.collapse.liver.dist1000.Robj", verbose=T)

N.sub.all <- subset(N.long, peak %in% S.collapse$peak)

# Each peak is a sample in 190 dimension space
# label each peak as "Liver" or "Flat"

peaks.hash <- hash(as.character(S.collapse$peak), as.character(S.collapse$peak.type))

# annotate N.sub.all
N.sub.all$peak.type <- sapply(as.character(N.sub.all$peak), function(p){
  return(peaks.hash[[p]])
})

N.mat <- dcast(data = N.sub.all, formula = peak ~ motif, 
               value.var = "sitecount", 
               fun.aggregate = sum)  # handle multple motifs in a peak by summing

N.mat$peak.type <- sapply(as.character(N.mat$peak), function(p){
  return(peaks.hash[[p]])
})

N.mat$peak <- NULL

fit <- lda(peak.type ~ ., data = N.mat)
plot(fit)
abline(a = 0, b = 1)
