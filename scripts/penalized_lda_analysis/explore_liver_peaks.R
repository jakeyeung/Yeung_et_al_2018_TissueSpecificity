# 2016-02-18
# After saving liver peaks, explore it to see what it looks like

library(penalizedLDA)
library(reshape2)
library(dplyr)
library(hash)

# Load --------------------------------------------------------------------

# load("Robjs/N.long.flatliver_genes.all_motifs.dist1000.Robj", verbose=T)
# load("Robjs/N.long.flatliver_genes.all_motifs.dist1000.Robj", verbose=T)
# load("Robjs/N.long.all_genes.all_signif_motifs.Robj", verbose=T)
load("Robjs/N.long.liverflat_genes.all_motifs.10000dist.Robj", verbose=T)

# load("Robjs/S.long.Robj", verbose=T)

# load("Robjs/S.liv.liverspec.peaks.Robj", verbose=T)  # fewer
load("Robjs/S.liver.liverspec.peaks.minamp02.Robj", verbose=T)  # more
S.sub.collapse.liv <- S.sub.collapse.filt
# load("Robjs/S.flat.liverspec.peaks.Robj", verbose=T)  # fewer
load("Robjs/S.flat.liverspec.peaks.top70percent.Robj", verbose=T)  # more
S.sub.collapse.flat <- S.sub.collapse.filt
rm(S.sub.collapse.filt)

# Filter peaks  -----------------------------------------------------------

# how many genes do we got?
(n.genes.liv <- length(unique(S.sub.collapse.liv$gene)))
(n.genes.flat <- length(unique(S.sub.collapse.flat$gene)))

peaks <- unique(c(as.character(S.sub.collapse.liv$peak), as.character(S.sub.collapse.flat$peak)))
N.long.filt.sub <- subset(N.long, peak %in% peaks)
str(N.long.filt.sub)


# Label model -------------------------------------------------------------

liv.genes <- unique(as.character(S.sub.collapse.liv$gene))
flat.genes <- unique(as.character(S.sub.collapse.flat$gene))

model.hash <- hash(c(liv.genes, flat.genes), c(rep("Liver", length(liv.genes)), rep("Flat", length(flat.genes))))

N.long.filt.sub$model <- sapply(as.character(N.long.filt.sub$gene), function(g) model.hash[[g]])

N.long.filt.sub$motif <- as.factor(N.long.filt.sub$motif)
N.long.filt.sub$peak <- as.factor(N.long.filt.sub$peak)
N.long.filt.sub$model <- as.factor(N.long.filt.sub$model)
N.long.filt.sub$sitecount <- as.numeric(N.long.filt.sub$sitecount)

# Run penalized LDA -------------------------------------------------------

# which motifs separate between Liv and Flat?
N.liv.mat <- dcast(data = subset(N.long.filt.sub, model == "Liver"), formula = peak ~ motif, value.var = "sitecount", fill = 0, fun.aggregate = sum)
N.flat.mat <- dcast(data = subset(N.long.filt.sub, model == "Flat"), formula = peak ~ motif, value.var = "sitecount", fill = 0, fun.aggregate = sum)

labs <- c(rep(1, nrow(N.liv.mat)), rep(2, nrow(N.flat.mat)))

common.motifs <- intersect(colnames(N.liv.mat), colnames(N.flat.mat))

N.livflat <- rbind(N.liv.mat[, common.motifs], N.flat.mat[, common.motifs])

rownames(N.livflat) <- N.livflat$peak; N.livflat$peak <- NULL

out <- PenalizedLDA(as.matrix(N.livflat), labs, lambda=0.1, K=1)

print(out)
plot(seq(length(out$discrim)), out$discrim)
text(seq(length(out$discrim)), out$discrim, , labels = colnames(N.livflat))


# Plot separation ---------------------------------------------------------

plot(out$xproj, out$y)
out.df <- data.frame(xproj = out$xproj, lab = out$y)
boxplot(xproj ~ lab, data = out.df)

# plot separation using a single dimension
jmotif <- "RORA.p2"
jmotif <- "ATF2.p2"
plot(labs, N.livflat[, jmotif])
