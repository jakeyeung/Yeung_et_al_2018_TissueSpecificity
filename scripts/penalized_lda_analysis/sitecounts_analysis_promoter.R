# 2016-03-03
# sitecounts_analysis_promoter.R
# Do on promoter (multiple background sets, summarize in x and y axis)

rm(list=ls())

library(ggplot2)
library(dplyr)
library(reshape2)
library(penalizedLDA)

# Function ----------------------------------------------------------------

source("scripts/functions/LdaFunctions.R")



# Load --------------------------------------------------------------------

load("Robjs/N.long.promoters_500.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)

# Stuff -------------------------------------------------------------------

# bg <- "Flat"
# fg.models <- "Adr"

# jgrep <- "^Kidney;Liver$|^Kidney,Liver$"
# fg.models <- unique(fits.best[grepl(jgrep, fits.best$model), ]$model)

fg.models <- c("Kidney;Liver", "Kidney,Liver")
# fg.models <- c("Liver")

fg.models <- c("Adr;Liver")

flat.models <- ""
rhyth.models <- as.character(subset(fits.best, n.rhyth >= 8)$model)

fg.genes <- subset(fits.best, model %in% fg.models)$gene
print(length(fg.genes))

fg.mat <- LongToMat(subset(N.long, gene %in% fg.genes & model %in% fg.models))
flat.mat <- LongToMat(subset(N.long, model %in% flat.models))
rhyth.mat <- LongToMat(subset(N.long, model %in% rhyth.models))

out.flat <- DoLdaPromoters(fg.mat, flat.mat)
out.rhyth <- DoLdaPromoters(fg.mat, rhyth.mat)


# Plot xy  ----------------------------------------------------------------

out.df <- data.frame(motif = colnames(out.flat$x), discrim.flat = out.flat$discrim, discrim.rhyth = out.rhyth$discrim)

ggplot(out.df, aes(x = discrim.flat, y = discrim.rhyth, label = motif)) + geom_text() + ggtitle(paste("Rhyth model:", paste(fg.models, collapse='|'))) + 
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
  theme_bw(24)

