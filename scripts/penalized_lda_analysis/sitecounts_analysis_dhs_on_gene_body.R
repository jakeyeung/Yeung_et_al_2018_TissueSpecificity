# 2016-02-19
# do penalized LDA using DHS on gene body

set.seed(0)

library(reshape2)
library(penalizedLDA)

# Load --------------------------------------------------------------------

load("Robjs/N.sitecounts.dhs.6.tissues.norm.Robj", verbose=T)  # N
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)

# Do LDA ------------------------------------------------------------------

minamp <- 0.25
liv.genes <- subset(fits.best, model == "Liver" & amp.avg > minamp)$gene

bkgrd.genes <- liv.genes  # same but take DHS from different tissues

N.liv <- subset(N, gene %in% liv.genes & tissue == "Liver")
N.bkgrd.all <- subset(N, gene %in% liv.genes & tissue != "Liver")

# sample to have same number as N.liv
(n.frgrd <- nrow(N.liv))
liv.genes.N <- unique(N.liv$gene)

# Run LDA -----------------------------------------------------------------

# make a new label gene;tissue which will be my rowname
N.liv$genetiss <- paste(N.liv$gene, N.liv$tissue, sep = ";")
N.bkgrd.all$genetiss <- paste(N.bkgrd.all$gene, N.bkgrd.all$tissue, sep = ";")

N.liv.mat <- dcast(data = N.liv, formula = genetiss ~ motif, fill = 0, value.var = "sitecount.norm")
N.bkgrd.all.mat <- dcast(N.bkgrd.all, genetiss ~ motif, fill = 0, value.var = "sitecount.norm")

labs <- c(rep(1, nrow(N.liv.mat)), rep(2, nrow(N.bkgrd.all.mat)))

N.mat.merged <- rbind(N.liv.mat, N.bkgrd.all.mat)
rownames(N.mat.merged) <- N.mat.merged$genetiss; N.mat.merged$genetiss <- NULL

out <- PenalizedLDA(as.matrix(N.mat.merged), labs, lambda = 0.1, K = 1, standardized = FALSE)

print(out)
plot(seq(length(out$discrim)), out$discrim)
text(seq(length(out$discrim)), out$discrim, , labels = colnames(N.mat.merged))

plot(out$xproj, out$y)
out.df <- data.frame(xproj = out$xproj, lab = out$y)
boxplot(xproj ~ lab, data = out.df)

# plot separation using a single dimension
jmotif <- "RXRG_dimer.p3"
jmotif <- "ONECUT1.2.p2"
out.df.motif <- data.frame(sitecount = N.mat.merged[, jmotif], lab = labs)
# plot(labs, N.mat.merged[, jmotif])
boxplot(sitecount ~ lab, data = out.df.motif)
