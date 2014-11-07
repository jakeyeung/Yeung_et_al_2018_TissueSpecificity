# Jake Yeung
# November 4 2014
# higher_order_svd.R
#
# Explore data using higher order single value decomposition
# use rTensor package

# Library and custom functions --------------------------------------------

library(rTensor)

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'hosvd.R'))  # modified from rTensor v 1.1 by James Li
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'PcaPlotFunctions.R'))  # for PCA and periodogram calcs

# Functions ---------------------------------------------------------------

ScaleTissue <- function(dat, N.TIMEPTS=24, N.TISSUES=12){
  # scale tissues (grouped by every N.TIMEPTS per tissue)
  # dat: expression matrix. Rows are genes. Columns are samples 
  # there are N.TIMEPTS + N.TISSUES columns.
  # grouped by tissues, ordered by time
  
#   # number of datapoints in matrix
#   N <- nrow(dat) * ncol(dat)
#   # create empty dataframe with same col and row names
#   scaled.dat <- as.data.frame(matrix(rep(NA, N), nrow=nrow(dat), ncol=ncol(dat)))
#   rownames(scaled.dat) <- rownames(dat)
#   colnames(scaled.dat) <- colnames(dat)
  for (i in 1:N.TISSUES){
    print(i)
    start.index <- (i - 1) * N.TIMEPTS + 1  # 1, 25, 49...
    end.index <- i * N.TIMEPTS  # 24, 48 ...
    dat[, start.index:end.index] <- t(scale(t(dat[, start.index:end.index])))
  }
  return(dat)
}

ScaleTime <- function(dat, N.TIMEPTS=24, N.TISSUES=12){
  # scale TIME
  # dat: expression matrix. Rows are genes. Columns are samples
  # there are N.TIMEPTS + N.TISSUES columns
  # grouped by tissues, so it's not as pretty to collect
  # all on same time point.
  # number of datapoints in matrix
  
#   N <- nrow(dat) * ncol(dat)
#   # create empty dataframe with same col and row names
#   scaled.dat <- as.data.frame(matrix(rep(NA, N), nrow=nrow(dat), ncol=ncol(dat)))
#   rownames(scaled.dat) <- rownames(dat)
#   colnames(scaled.dat) <- colnames(dat)
  
  for (i in 1:N.TIMEPTS){
    print(i)
    indices <- seq(i, N.TIMEPTS * N.TISSUES, N.TIMEPTS)
    dat[, indices] <- t(scale(t(dat[, indices])))
  }
  return(dat)
}


# Test functions work -----------------------------------------------------

# dat.test <- as.data.frame(matrix(1:288*5, ncol=288, nrow=5))
# dat.test.tiss <- ScaleTissue(dat.test, 24, 12)
# print(dat.test[1, 1:24])
# print(dat.test.tiss[1, 1:24])

# dat.test <- as.data.frame(matrix(1:288*5, ncol=288, nrow=5))
# dat.test.time <- ScaleTime(dat.test)
# dat.test.time.tissue <- ScaleTissue(dat.test.time)
# 
# i <- seq(1, 288, 24)
# print(rowMeans(dat.test.time.tissue[1, 1:24]))
# print(rowMeans(dat.test.time.tissue[1, i]))

# Define constants --------------------------------------------------------

N.TIMEPTS <- 24  # 24 time points, 2 hours per time point (over 48 hrs)
N.TISSUES <- 12  # 12 tissues

# Load data ---------------------------------------------------------------


# define dirs
data_dir <- "data"
fname <- "hogenesch_2014_rma.txt"    # data reprocessed by RMA package
# fname <- "hogenesch_2014_rma.ensemblnames.txt"

# load data
data_path <- file.path(data_dir, fname)
print(paste("Reading data from,", data_path, "May take a few a minutes."))
dat <- read.table(data_path)
print("Read data to memory.")


# Get Colnames ------------------------------------------------------------

colnames(dat) <- ShortenSampNames(colnames(dat), show="tissue.time")
dat.colnames <- colnames(dat)  # in case I lose it later


# Scale by tissue and time ------------------------------------------------

# if you're looking at tissue, you scale by time.
dat.tissue <- ScaleTime(dat, N.TIMEPTS, N.TISSUES)

# if you're looking at time, you scale by tissue
dat.time <- ScaleTissue(dat, N.TIMEPTS, N.TISSUES)


# Convert to Tensor -------------------------------------------------------

# 3-mode tensor.
# mode 1: genes (N ~ tens of thousands)
# mode 2: number of time points (N = 24)
# mode 3: number of tissues (N = 12)

tnsr.tissue <- fold(as.matrix(dat.tissue), rs=1, cs=c(2, 3), modes=c(nrow(dat), N.TIMEPTS, N.TISSUES))
tnsr.time <- fold(as.matrix(dat.time), rs=1, cs=c(2, 3), modes=c(nrow(dat), N.TIMEPTS, N.TISSUES))

# sanity check
# Below is matrix exprs of row 1 gene. Rows are time. Columns are tissues.
# tnsr.tissue[1, 1:N.TIMEPTS, 1:N.TISSUES]


# Higher Order SVD --------------------------------------------------------

# HOSVD is not unique
tnsr.tissue.hosvd <- hosvd(tnsr.tissue, rank=c(N.TIMEPTS * N.TISSUES, N.TIMEPTS, N.TISSUES))  # ~few minutes
tnsr.time.hosvd <- hosvd(tnsr.time, rank=c(N.TIMEPTS * N.TISSUES, N.TIMEPTS, N.TISSUES))  # few minutes

# CP is often unique, but columns of U matrices are not orthnormal
# tnsr.tissue.cp <- cp(tnsr.tissue, num_components=12)  # WARNING: kind of slow
# tnsr.time.cp <- cp(tnsr.time, num_components=24)  # WARNING: SLOW +30 MIN

# print dimensions of each orthogonal matrix, U
lapply(tnsr.tissue.hosvd$U, dim)
lapply(tnsr.time.hosvd$U, dim)


# Plot PCAs ---------------------------------------------------------------

# tnsr.hosvd$U contains a list of matrices. 
# tnsr.hosvd$U[[1]] is size 35556 by 288  # "eigensamples?"
# tnsr.hosvd$U[[2]] is size 24 by 24  # "time-component eigengenes?"
# tnsr.hosvd$U[[3]] is size 12 by 12  # "tissue-component eigengenes?"

# take eigengenes from HOSVD 
eigengenes.time <- tnsr.time.hosvd$U[[2]]
eigengenes.tissue <- tnsr.tissue.hosvd$U[[3]]
# take eigengenes from CP
# eigengenes.time <- tnsr.time.cp$U[[2]]
# eigengenes.tissue <- tnsr.tissue.cp$U[[3]]

dim(eigengenes.time)
dim(eigengenes.tissue)


# Plot tissue PCAs 2D ------------------------------------------------------


for (i in 1:10){
  x <- eigengenes.tissue[, i]
  y <- eigengenes.tissue[, i + 1]
  plot(x, y, main=paste0('TISSUE: PCA ', i, ' vs. ', 'PCA ', i + 1))
  labels <- dat.colnames[seq(1, 288, 24)]
  text(x, y, labels, pos=3)
}




# Plot time PCAs 2D --------------------------------------------------------


# colors for time. Repeat twice for 24 hour oscillations over 48 hrs
cols <- rep(c('black', 'red', 'cyan', 'blue', 'lawngreen', 'burlywood'), 2)
for (i in 1:10){
  x <- eigengenes.time[, 1]
  y <- eigengenes.time[, i + 1]
  plot(x, y, main=paste0('TIME: PCA ', 1, ' vs. ', 'PCA ', i + 1), col=cols)
  labels <- rep(seq(0, 22, 2), 2)
  text(x, y, labels, pos=3, col=cols)
}


# Plot 1D PCA plus periodogram --------------------------------------------


# We will fit one linear models:
# 1) linear model with common intercept and common oscillating component
# and common oscillating component

N <- 24  # number of samples == number of time periods in this case.
T <- 24  # 24 hours in a period

# Fit PCA component of interest: loop to try many different PCAs
# user changeable range
for (pca_vector in 1:10) {
  # Create response vector, which is loadings
  
  y <- eigengenes.time[, pca_vector]
  
  # BEGIN: plot periodograms to see which frequency has high activity
  freq.and.periodogram <- CalculatePeriodogram(y)  # returns a list
  freq <- freq.and.periodogram$freq
  periodogram <- freq.and.periodogram$p.scaled
  periodogram.unscaled <- freq.and.periodogram$p.unscaled
  
  # Calculate top 5 frequencies
  max.freqs <- FindMaxFreqs(freq, periodogram)
  
  PlotPeriodogram(freq, periodogram, title=paste("Periodogram for PCA component:", pca_vector))
  # add vertical line at max frequency
  max.f <- max.freqs[1]
  # calculate period from frequency
  max.T <- (1 / max.f) * 2  # multiply by 2 because samples are every 2 hours 
  abline(v=max.f, col='blue', lwd=2)
  # add text to show freq and period.
  # x offset of +0.02 so you can see the text
  text(max.f + 0.02, 0, paste0("T=", signif(max.T, digits=3), "hrs"))
  
  PlotLoadings(y, title=paste("Vector Loadings for PCA component:", pca_vector))
  # END: plot periodograms to see which frequency has high activity
  
  # BEGIN: Linear fit for period of 24 hours
  # Create sequence of t = [0, 2, 4, ... (N - 1)]
  t <- seq(0, 2 * N - 1, 2)
  
  # set my angular frequency
  omega <- (2 * pi) / (T)
  
  # fit my lm using cos and sin with angular frequency
  fit <- lm(y ~ cos(omega * t) + sin(omega * t))
  # END: Linear fit for period of 24 hours
  
  # Print some statementsy to describe fit
  cat("*********************************\n")
  cat(paste0("PCA ", pca_vector, "\n"))
  cat("\n")
  
  cat("Max freqs\n")
  print(max.freqs)
  
  cat("Simple fit:\n")
  print(anova(fit))
  cat("\n")
  
}




