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
source("scripts/functions/RemoveP2Name.R")


# Load --------------------------------------------------------------------

load("Robjs/N.long.promoters_500.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)

# Stuff -------------------------------------------------------------------

# bg <- "Flat"
fg.models <- "Adr"

# jgrep <- "^Kidney;Liver$|^Kidney,Liver$"
# fg.models <- unique(fits.best[grepl(jgrep, fits.best$model), ]$model)

# fg.models <- c("Kidney;Liver", "Kidney,Liver")
# fg.models <- c("Liver")

# fg.models <- c("Adr;Liver")

flat.models <- ""
rhyth.models <- as.character(subset(fits.best, n.rhyth >= 8)$model)

fg.genes <- subset(fits.best, model %in% fg.models)$gene
print(length(fg.genes))

fg.mat <- LongToMat(subset(N.long, gene %in% fg.genes & model %in% fg.models))
flat.mat <- LongToMat(subset(N.long, model %in% flat.models))
rhyth.mat <- LongToMat(subset(N.long, model %in% rhyth.models))

# SINGLE FACTORS
mat.fgbg <- bind_rows(fg.mat, flat.mat)
mat.fgbg[is.na(mat.fgbg)] <- 0
rownames(mat.fgbg) <- c(rownames(fg.mat), rownames(flat.mat))
labels <- c(rep(1, nrow(fg.mat)), rep(2, nrow(flat.mat)))
# remove columns with 0 variance 
mat.fgbg[which(colSums(mat.fgbg) == 0)] <- list(NULL)
out <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  
PlotLdaOut(out)

# CROSS PRODUCTS
# do crosses
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)
mat.fgbg.cross <- CrossProduct(mat.fgbg, remove.duplicates = TRUE)
dim(mat.fgbg.cross)
mat.fgbg.cross <- cbind(mat.fgbg, mat.fgbg.cross)
# remove columns with 0 variance 
mat.fgbg.cross[which(colSums(mat.fgbg.cross) == 0)] <- list(NULL)
out.cross <- PenalizedLDA(mat.fgbg.cross, labels, lambda = 0.01, K = 1, standardized = FALSE)
m <- SortLda(out.cross)
print(length(m))
BoxplotLdaOut(out.cross)
PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE)

out.flat <- DoLdaPromoters(fg.mat, flat.mat)
out.rhyth <- DoLdaPromoters(fg.mat, rhyth.mat)


# Plot xy  ----------------------------------------------------------------

out.df <- data.frame(motif = colnames(out.flat$x), discrim.flat = out.flat$discrim, discrim.rhyth = out.rhyth$discrim)

ggplot(out.df, aes(x = discrim.flat, y = discrim.rhyth, label = motif)) + geom_text() + ggtitle(paste("Rhyth model:", paste(fg.models, collapse='|'))) + 
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
  theme_bw(24)

