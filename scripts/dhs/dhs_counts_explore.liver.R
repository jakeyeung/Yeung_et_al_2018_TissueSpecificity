# Look at DHS counts data
# Date: 2015-03-17
rm(list=ls())

library(ggplot2)

# Functions ---------------------------------------------------------------

GetReplicates <- function(cnames, get.index=TRUE){
  # from colnames of data, retrieve replicate names
  # expects that colnames are not called "chr", "start", "end", "zscore" or "total"
  # get.index -> return the index from cnames otherwise get names
  colnames.remove <- c("chr", "start", "end", "zscore", "total")
  if (get.index){
    return(which(!cnames %in% colnames.remove))
  }
  else {
    return(cnames[which(!cnames %in% colnames.remove)])
  }
}

fit <- function(x, a, b){
  return(exp(a*x + b))
}

# Main --------------------------------------------------------------------


# Get data ----------------------------------------------------------------


dhs.path <- "~/projects/tissue-specificity/data/beds/ENCODE:Seq:Liver/combined.nogc.bedfixed"
# dhs.path <- "/home/yeung/projects/tissue-specificity/data/beds/ENCODE:Seq:Kidney/combined.bed"
dhs.dat <- read.table(dhs.path, header = TRUE, sep = "\t")  # gigantic file ~ 1.2 gigs in memory
save(dhs.dat, file = "~/projects/tissue-specificity/docs/2015-03-19-felix/dhs_dat.Robj")
jtissue <- "Liver"

# subset constants --------------------------------------------------------
set.seed(0)
n.sub <- 0.01 * nrow(dhs.dat)
rows.sub <- sample(seq(1:nrow(dhs.dat)), n.sub)
i.reps <- GetReplicates(colnames(dhs.dat))
n.samps <- length(i.reps)


# Normalize by total counts -----------------------------------------------

dhs.reps <- dhs.dat[, i.reps]
sum.samps <- apply(dhs.reps, 2, sum)
dhs.reps <- dhs.reps / sum.samps * 10^6

barplot(sum.samps, las = 2)

# Plot densities  ---------------------------------------------------------

dhs.reps.liversub <- dhs.reps[rows.sub, ]
for (i in 1:ncol(dhs.reps.liversub)){
  plot(density(log2(dhs.reps[rows.sub, i])), main = paste("rep", i))
}

# Identify outlier --------------------------------------------------------

samp.outlier <- "UwStam_mLiver.DS16858.FC62FJY.4_004"
samp.normal <- "wgEncodeUwDnaseLiverC57bl6MAdult8wksRawDataRep1"
samp.normal2 <- "wgEncodeUwDnaseLiverC57bl6MAdult8wksRawDataRep2"

# plot pairs
n.sub.forpairs <- 0.001 * nrow(dhs.reps)
rows.subpair <- sample(seq(1:nrow(dhs.reps)), n.sub.forpairs)
# this takes super long and prone to crashes
# pdf("plots/dhs/pair_plots_scattermatrix.pdf")
# pairs(log2(dhs.dat[rows.subpair, i.reps]))  # UwStam_mLiver.DS16858.FC62FJY.4_004 is an outlier
# dev.off()

# Plot pairwise between outlier and normal and normal-normal
pairs(log2(dhs.reps[rows.subpair, c(samp.outlier, samp.normal), ]))
pairs(log2(dhs.reps[rows.subpair, c(samp.normal2, samp.normal), ]))
pairs((log2(dhs.reps[rows.subpair, c("UwStam.mLiver.DS16740.FC62FJL.7", "UwStam_mLiver.DS16853.FC62FJY.3_004"), ])))

# Plot density distributions for outlier and two normals
counts.outlier <- dhs.reps[, samp.outlier][which(dhs.reps[, samp.outlier] > 0)]
counts.normal <- dhs.reps[, samp.normal][which(dhs.reps[, samp.normal] > 0)]
counts.normal2 <- dhs.reps[, samp.normal2][which(dhs.reps[, samp.normal2] > 0)]
qplot(x = counts.outlier, geom = "density") + scale_x_log10() + ggtitle(paste0("Sample ", samp.outlier))
qplot(x = counts.normal, geom = "density") + scale_x_log10() + ggtitle(paste0("Sample ", samp.normal))
qplot(x = counts.normal2, geom = "density") + scale_x_log10() + ggtitle(paste0("Sample ", samp.normal2))

# Remove outlier ----------------------------------------------------------

dhs.dat.no.out <- dhs.reps
dhs.dat.no.out[, samp.outlier] <- NULL

# get mean and sd ---------------------------------------------------------

epsilon <- 10^-3
sums.reps <- apply(dhs.dat.no.out[rows.sub, ], 1, sum)
# dhs.dat.no.out <- scale(dhs.dat.no.out)
means.reps <- apply(dhs.dat.no.out[rows.sub, ], 1, mean)
sds.reps <- apply(dhs.dat.no.out[rows.sub, ], 1, sd)

qplot(x = sums.reps + 1, geom = "density") + scale_x_log10()

ggplot(data = data.frame(Mean = means.reps, SD = sds.reps), 
       aes(x = Mean, y = SD)) + 
  geom_point(alpha = 0.05) + 
  scale_x_log10() + 
  scale_y_log10() +
  ggtitle("Standard deviation vs mean")

ggplot(data = data.frame(Mean = means.reps, CV = sds.reps / means.reps), 
       aes(x = Mean, y = CV)) + 
  geom_point(alpha = 0.05) + 
  scale_x_log10() + 
  scale_y_log10() +
  ggtitle("Coefficient of variation vs mean")


# # Fit functions -----------------------------------------------------------
# 
# # fitting code
# y = sds.log2.reps / means.log2.reps
# x = means.log2.reps
# 
# fitModel = nls(y ~ fit(x, a, b), start = list(a = -1, b = 1))
# params = coef(fitModel)
# 
# plot(log2(means.reps + 1), log2(sds.reps + 1), main = paste0("SD vs Mean: ", jtissue), xlab = "log2(Mean)", ylab = "log2(SD)")
# 
# # Plot log CV against log counts
# plot(log2(means.reps / sds.reps + 1), log2(means.reps), main = paste0("CV vs mean: ", jtissue), xlab = "log2(mean)", ylab = "log2(CV)")

