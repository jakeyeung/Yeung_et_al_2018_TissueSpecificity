# 2016-02-19
# do penalized LDA using DHS on gene body
# redo using a new set of sitecounts matrix derived from get_tissue_spec_peaks.R

rm(list=ls())

library(reshape2)
library(penalizedLDA)

# Load --------------------------------------------------------------------

use.cross <- TRUE

if (!use.cross){
  load("Robjs/N.long.sum.bytiss.Robj", verbose=T)
} else {
  load("Robjs/N.long.sum.bytiss.RXRGcross.Robj", verbose=T)
}

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)

# Do LDA ------------------------------------------------------------------

minamp <- 0.25
liv.genes <- subset(fits.best, model == "Liver" & amp.avg > minamp)$gene

bkgrd.genes <- liv.genes  # same but take DHS from different tissues

all.tiss <- c('Cere', 'Heart', 'Kidney', 'Liver', 'Lung', 'Mus')

# fg.tiss <- c("Liver", "Kidney")
# fg.tiss <- c("Liver")
# bg.tiss <- all.tiss[which(! all.tiss %in% fg.tiss)]

fg.tiss <- c("Liver")
bg.tiss <- all.tiss[which(! all.tiss %in% fg.tiss & all.tiss != "Kidney")]

# fg.tiss <- c("Liver")
# bg.tiss <- c("Kidney")

N.liv <- subset(N.long.sum.bytiss, gene %in% liv.genes & tissue %in% fg.tiss)
N.bkgrd.all <- subset(N.long.sum.bytiss, gene %in% liv.genes & tissue %in% bg.tiss)

# sample to have same number as N.liv
(n.frgrd <- nrow(N.liv))
liv.genes.N <- unique(N.liv$gene)

# Run LDA -----------------------------------------------------------------

# make a new label gene;tissue which will be my rowname
N.liv$genetiss <- paste(N.liv$gene, N.liv$tissue, sep = ";")
N.bkgrd.all$genetiss <- paste(N.bkgrd.all$gene, N.bkgrd.all$tissue, sep = ";")

if (!use.cross){
  sitecount.name <- "sitecount.max"
} else {
  sitecount.name <- "sitecount.cross.max"
}

N.liv.mat <- dcast(data = N.liv, formula = genetiss ~ motif, fill = 0, value.var = sitecount.name)
N.bkgrd.all.mat <- dcast(N.bkgrd.all, genetiss ~ motif, fill = 0, value.var = sitecount.name)

labs <- c(rep(1, nrow(N.liv.mat)), rep(2, nrow(N.bkgrd.all.mat)))

N.mat.merged <- rbind(N.liv.mat, N.bkgrd.all.mat)
rownames(N.mat.merged) <- N.mat.merged$genetiss; N.mat.merged$genetiss <- NULL

out <- PenalizedLDA(as.matrix(N.mat.merged), labs, lambda = 0.1, K = 1, standardized = FALSE)

print(out)
discrim.filt <- sort(out$discrim[which(out$discrim != 0)], decreasing = FALSE)
cnames <- colnames(N.mat.merged)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
plot(seq(length(discrim.filt)), discrim.filt)
text(seq(length(discrim.filt)), discrim.filt, labels = cnames)

print(cnames[1:20])

# 
# plot(out$xproj, out$y)
# out.df <- data.frame(xproj = out$xproj, lab = out$y)
# boxplot(xproj ~ lab, data = out.df)
# 
# # plot separation using a single dimension
# jmotif <- "RXRG_dimer.p3"
# jmotif <- "ONECUT1,2.p2"
# jmotif <- "HNF1A.p2"
# jmotif <- "FOXA2.p3"
# jmotif <- "MAFB.p2"
# 
# out.df.motif <- data.frame(sitecount = N.mat.merged[, jmotif], lab = labs)
# boxplot(sitecount ~ lab, data = out.df.motif, main = jmotif)
# # plot(labs, N.mat.merged[, jmotif])


# Where are the RXRG_dimer motifs located? --------------------------------


