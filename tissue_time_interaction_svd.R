# Jake Yeung
# tissue_time_interactions_svd.R
# 2014-11-07
# Fix periodicity in time component. SVD on tissues.


# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'PcaPlotFunctions.R'))  # for PCA and periodogram calcs
source(file.path(functions.dir, 'FourierFunctions.R'))  # for Fourier stuff

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

# RHYTHMIC
OMEGA <- 2 * pi / PERIOD
dat.time.projected <- ProjectToPeriodicTime(as.matrix(dat), 
                                            N.TISSUES, 
                                            N.TIMEPTS, 
                                            INTERVAL, 
                                            OMEGA, 
                                            dat.tissuenames)

Perform PCA -------------------------------------------------------------

# We'll do this manually with SVD decomposition.
dat.pca <- svd(scale(t(dat.time.projected)))
plot(dat.pca$d^2)  # Manual screeplot. 
print(lapply(dat.pca, dim))  # $u is 12x12, $v is k-genes x 12

# You can check that this is equivalent to prcomp
# dat.pca.check <- prcomp(scale(t(dat.time.projected)))
# screeplot(dat.pca.check, npcs=length(dat.pca.check$sdev), type="lines")


# Name U and V rows with samples and genes -------------------------------

rownames(dat.pca$u) <- dat.tissuenames
rownames(dat.pca$v) <- rownames(dat)

# Plot PCA ----------------------------------------------------------------

for (i in 1:5){
  # plotting real values: PCA 1 vs PCA 2 show Liver is different.
  # dat.pca$u has no Imaginary part? Arg are all either pi or 0...
  x <- dat.pca$u[, i]
  y <- dat.pca$u[, i + 1]
  plot(x, y, main=paste0('TISSUE: PCA ', i, ' vs. ', 'PCA ', i + 1))
  labels <- dat.tissuenames
  text(x, y, labels, pos=3)
}

# Find "modules" of important genes ---------------------------------------

# limit to top N genes that give highest loadings...
N.genes <- 500

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
        scale="none",
        main=paste('Original data'))








