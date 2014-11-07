# Jake Yeung
# tissue_time_interactions_svd.R
# 2014-11-07
# Fix periodicity in time component. SVD on tissues.


# Functions ---------------------------------------------------------------

source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'PcaPlotFunctions.R'))  # for PCA and periodogram calcs

ProjectToPeriodicTime <- function(Y , N.TISSUES, N.TIMEPTS, INTERVAL, OMEGA, col.names){
  # Project matrix of times to frequency domain
  # ARGS:
  #   Y: matrix of expression. Rows are genes. Columns contain tissues and time.
  #      In this case, expect tissues clustered together, ordered by time.
  #   col.names: column names of output Y's projected onto time. FALSE means no col.names
  #   N.TISSUES: number of tissues
  #   N.TIMEPTS: number of time points per tissue
  #   INTERVAL: interval between time points. e.g. 2 if sampled every 2hrs
  #   OMEGA: 2 * pi / PERIOD. Angular frequency. If omega = 0, matrix Y is 
  #   normalized by time components (not oscillating)
  # 
  # RETURNS:
  #   Y.time.projected: matrix of expression, projected onto the temporal axis.
  
  # track number of genes for dimension purposes
  N.GENES <- nrow(Y)
  # get times vector
  times.vec <- seq(length.out = N.TIMEPTS, by = INTERVAL)
  
  # init output matrix
  Y.time.projected <- matrix(NA, nrow=N.GENES, ncol=N.TISSUES)
  
  # identify row and colnames
  rownames(Y.time.projected) <- rownames(Y)  # same row names
  colnames(Y.time.projected) <- col.names  # from args  
  
  # BEGIN: project onto temporal axis
  for (i in 1:N.TISSUES){
    # get tissue.i across time
    index.start <- (i - 1) * N.TIMEPTS + 1  # 1, 25, 49...
    index.end <- i * N.TIMEPTS
    Y.tissue.i <- Y[, index.start:index.end]  # all genes
    # project tissues onto temporal axis
    # if OMEGA = 0, it is equivalent to getting average of Y.tissue.i
    T <- ((exp(1i * OMEGA * times.vec)) / (N.TIMEPTS))
    Y.time.projected[, i] <- Y.tissue.i %*% (exp(1i * OMEGA * times.vec)) / (N.TIMEPTS)
  }
  return(Y.time.projected)
}

GetTopNValues <- function(x, N){
  # Return top N values from vector x
  # 
  # ARGS:
  # x: vector
  # N: top N values to return
  # 
  # RETURNS:
  # x.top: list with components
  #   vals: top N values
  #   i: indices of the top N values
  
  # init list
  x.top <- list(vals=NULL, i=NULL)
  
  x.top$vals <- head(sort(x, decreasing=TRUE), N)
  x.top$i <- which(x %in% x.top$vals)  # not ordered!
  return(x.top)
}

# Test function works -----------------------------------------------------

# test ProjectToPeriodicTIme
#
# ROWS <- 3
# COLS <- 9
# TIMEPTS <- 3  # 3 time points per tissue
# TISSUES <- COLS / TIMEPTS
# Y <- matrix(seq(3 * 9), nrow=ROWS, ncol=COLS)
# out.colnames <- make.unique(rep('COL', TISSUES))
# rownames(Y) <- make.unique(rep('ROW', ROWS))
# (Y.t <- ProjectToPeriodicTime(Y, N.TISSUES=TISSUES, N.TIMEPTS=TIMEPTS, INTERVAL=3, OMEGA=0, out.colnames))

# # test GetTopNValues
# x <- seq(from=50, to=1)
# x.top <- GetTopNValues(x, 5)

# Define constants --------------------------------------------------------

N.TIMEPTS <- 24  # 24 time points, 2 hours per time point (over 48 hrs)
N.TISSUES <- 12  # 12 tissues
PERIOD <- 24  # 24 hours
INTERVAL <- 2  # hours between samples


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
dat.tissuenames <- dat.colnames[seq(1, 288, 24)]
# tissue names are Adr18... remove the 18 from tissue names
dat.tissuenames <- unname(sapply(dat.tissuenames, function(x){
  return(substr(x, 1, nchar(x) - 2))
}))


# Project onto flat and rhythmic time components --------------------------

OMEGA <- 2 * pi / PERIOD
dat.time.projected <- ProjectToPeriodicTime(as.matrix(dat), 
                                            N.TISSUES, 
                                            N.TIMEPTS, 
                                            INTERVAL, 
                                            OMEGA, 
                                            dat.tissuenames)



# Perform PCA -------------------------------------------------------------

# check SVD and prcomp are spitting out same matrices
# dat.pca.check <- prcomp(scale(t(dat.time.projected)))
# screeplot(dat.pca.check, npcs=length(dat.pca.check$sdev), type="lines")

dat.pca <- svd(scale(t(dat.time.projected)))
plot(dat.pca$d^2)  # screeplot
lapply(dat.pca, dim)


# Name U and V rows with samples and genes -------------------------------

rownames(dat.pca$u) <- dat.tissuenames
rownames(dat.pca$v) <- rownames(dat)

# Plot PCA ----------------------------------------------------------------

for (i in 1:5){
  x <- dat.pca$u[, i]
  y <- dat.pca$u[, i + 1]
  plot(x, y, main=paste0('TISSUE: PCA ', i, ' vs. ', 'PCA ', i + 1))
  labels <- dat.tissuenames
  text(x, y, labels, pos=3)
}


# Find "modules" of important genes ---------------------------------------

# limit to top N genes that give highest loadings...
N.genes <- 200

# get top N values from V. Therefore filter...
# Value calculated from Modulus of elements in vector to take into account
# both Re and Im

for (PCA in 1:10){
  top.N <- GetTopNValues(Mod(dat.pca$v[, PCA]), N.genes)
  
  # filter expression matrix for only top contributing genes
  # recreate approximation of gene exprs with PCA
  v.truncated <- dat.pca$v[top.N$i, ]  # take v matrix, truncate for heatmap
  d.row <- rep(0, length(dat.pca$d))
  d.row[PCA] <- dat.pca$d[PCA]  # create row vector, all zeros except at diagonal
  
  # approximate gene exprs with PCA. Multiply first diagonal with v, then multiply that by u
  d.v <- d.row %*% t(v.truncated)
  dat.by.pca <- dat.pca$u[, PCA] %*% d.v
  # transpose it so genes are rows, tissues are columns
  dat.by.pca <- t(dat.by.pca)  # N.genes by N.tissues matrix
  
  colnames(dat.by.pca) <- dat.tissuenames
  
  # plot dat heatplot in Amplitude and Phase
  # phase
  heatmap(as.matrix(Arg(dat.by.pca)), 
          Rowv=NA, 
          Colv=NA, 
          main=paste('PCA:', PCA))  
}
# our original gene expression...
dat.filtered <- dat.time.projected[top.N$i, ]
heatmap(as.matrix(Arg(dat.filtered)), 
        Rowv=NA, 
        Colv=NA,
        main=paste('Original data'))








