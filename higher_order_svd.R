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

# Functions ---------------------------------------------------------------

ScaleTissues <- function(dat, N.TIMEPTS=24, N.TISSUES=12){
  # scale tissues (grouped by every N.TIMEPTS per tissue)
  # dat: expression matrix. Rows are genes. Columns are samples 
  # there are N.TIMEPTS + N.TISSUES columns.
  # grouped by tissues, ordered by time
  
  # number of datapoints in matrix
  N <- nrow(dat) * ncol(dat)
  # create empty dataframe with same col and row names
  scaled.dat <- as.data.frame(matrix(rep(NA, N), nrow=nrow(dat), ncol=ncol(dat)))
  rownames(scaled.dat) <- rownames(dat)
  colnames(scaled.dat) <- colnames(dat)
  
  for (row.i in 1:nrow(dat)){
    row.vec <- dat[row.i, ]
    for (i in 1:N.TISSUES){
      start.index <- (i - 1) * N.TIMEPTS + 1  # start at 1, 25, 49 ...
      end.index <- i * N.TIMEPTS  # ends at 24, 48, ...
      # get vector of 1 tissue across all times
      tissue.vec <- row.vec[start.index:end.index]
      # scale, mean centre and scale
      tissue.vec.scaled <- as.vector(scale(t(tissue.vec), scale=FALSE))
      # update row.vec with new scaled value
      row.vec[start.index:end.index] <- tissue.vec.scaled
    }
    scaled.dat[row.i, ] <- row.vec
  }
  return(scaled.dat)
}

# test function
dat.test <- as.data.frame(matrix(1:288*5, nrow=5, ncol=288))
dat.test.tiss <- ScaleTissues(dat.test, 24, 12)



ScaleTime <- function(dat, N.TIMEPTS, N.TISSUES){
  # scale TIME
  # dat: expression matrix. Rows are genes. Columns are samples
  # there are N.TIMEPTS + N.TISSUES columns
  # grouped by tissues, so it's not as pretty to collect
  # all on same time point.
  scaled.dat <- apply(dat, 1, function(row.vec){
    for (i in 1:N.TIMEPTS){
      # if grouped by tissues, to get time, we need
      # to grab element at every 24 elements
      indices <- seq(i, N.TIMEPTS * N.TISSUES, N.TIMEPTS)
      time.vec <- row.vec[indices]
      # scale
      time.vec.scaled <- scale(time.vec)
      # update row.vec
      row.vec[indices] <- time.vec.scaled
    }
    return(row.vec)
  })
  str(scaled.dat)
  return(scaled.dat)
}


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

colnames(dat) <- ShortenSampNames(colnames(dat), show="tissue")
dat.colnames <- colnames(dat)  # in case I lose it later


# Scale by tissue and time ------------------------------------------------

print(row.vec <- dat[1, ])

dat <- ScaleTissues(dat, N.TIMEPTS, N.TISSUES)

print(row.vec.tissue.scaled <- dat[1, ])

dat <- ScaleTime(dat, N.TIMEPTS, N.TISSUES)

print(row.vec.time.scaled <- dat[1, ])

# Convert to Tensor -------------------------------------------------------

# 3-mode tensor.
# mode 1: genes (N ~ tens of thousands)
# mode 2: number of time points (N = 24)
# mode 3: number of tissues (N = 12)

tnsr <- fold(dat, rs=1, cs=c(2, 3), modes=c(nrow(dat), N.TIMEPTS, N.TISSUES))

# sanity check
# Below is matrix exprs of row 1 gene. Rows are time. Columns are tissues.
# tnsr[1, 1:N.TIMEPTS, 1:N.TISSUES]


# Higher Order SVD --------------------------------------------------------

# HOSVD, truncate after 100 because I get subscript out of bound errors?!

tnsr.hosvd <- hosvd(tnsr, rank=c(N.TIMEPTS * N.TISSUES, N.TIMEPTS, N.TISSUES))

# print dimensions of each orthogonal matrix, U
lapply(tnsr.hosvd$U, dim)


# Plot PCAs ---------------------------------------------------------------

# tnsr.hosvd$U contains a list of matrices. 
# tnsr.hosvd$U[[1]] is size 35556 by 288  # "eigensamples?"
# tnsr.hosvd$U[[2]] is size 24 by 24  # "time-component eigengenes?"
# tnsr.hosvd$U[[3]] is size 12 by 12  # "tissue-component eigengenes?"

eigengenes.time <- tnsr.hosvd$U[[2]]
eigengenes.tissue <- tnsr.hosvd$U[[3]]

plot(eigengenes.time[, 1], eigengenes.time[, 2])
plot(eigengenes.tissue[, 1], eigengenes.tissue[, 2])


