# Jake Yeung
# tissue_time_interactions_svd.R
# 2014-11-07
# Fix periodicity in time component. SVD on tissues.


# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'PcaPlotFunctions.R'))  # for PCA and periodogram calcs
source(file.path(functions.dir, 'FourierFunctions.R'))  # for Fourier stuff
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data

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
fname <- "hogenesch_2014_rma.genenames.txt"  # has duplicate gene names

# load data
data_path <- file.path(data_dir, fname)
print(paste("Reading data from,", data_path, "May take a few a minutes."))
dat <- read.table(data_path, header=TRUE)
print("Read data to memory.")



# If needed handle duplicate row names ------------------------------------

rownames(dat) <- make.names(dat$gene, unique=TRUE)
# genes <- array.exprs$gene
# coords <- array.exprs$coordinates
# probeid <- array.exprs$ID_REF

drop.cols <- c("gene")
dat <- dat[, !(names(dat) %in% drop.cols)]

Peek(dat)  # expect gene names as row names, tissues in columns

# Get Colnames ------------------------------------------------------------

colnames(dat) <- ShortenSampNames(colnames(dat), show="tissue.time")
dat.colnames <- colnames(dat)  # in case I need it later
dat.tissuenames <- dat.colnames[seq(1, 288, 24)]
# tissue names are Adr18... remove the 18 from tissue names
dat.tissuenames <- unname(sapply(dat.tissuenames, function(x){
  return(substr(x, 1, nchar(x) - 2))
}))

Peek(dat)  # sample names are more readable...

# Project onto flat and rhythmic time components --------------------------

# RHYTHMIC
OMEGA <- 2 * pi / PERIOD
dat.time.projected <- ProjectToPeriodicTime(as.matrix(dat), 
                                            N.TISSUES, 
                                            N.TIMEPTS, 
                                            INTERVAL, 
                                            OMEGA, 
                                            dat.tissuenames)



# PCA ---------------------------------------------------------------------

# SHOULD WE SCALE OR NOT? 
# dat.pca <- svd(t(scale(t(dat.time.projected))))  # scaled
dat.pca <- svd(dat.time.projected)  # not scaled
plot(dat.pca$d^2, type='o')  # Manual screeplot. 
print(lapply(dat.pca, dim))  # $u is 12x12, $v is k-genes x 12

# You can check that this is equivalent to prcomp
# dat.pca.check <- prcomp(scale(t(dat.time.projected)))
# screeplot(dat.pca.check, npcs=length(dat.pca.check$sdev), type="lines")
# summary(dat.pca.check)

# Name U and V rows with samples and genes -------------------------------

rownames(dat.pca$v) <- dat.tissuenames
rownames(dat.pca$u) <- rownames(dat)

# Plot PCA ----------------------------------------------------------------

for (i in 1:10){
  x <- Arg(dat.pca$v[, i])
  y <- Arg(dat.pca$v[, i + 1])
  plot(x, y, main=paste0('TISSUE: PCA ', i, ' vs. ', 'PCA ', i + 1))
  labels <- dat.tissuenames
  text(x, y, labels, pos=3)
}

# Find "modules" of important genes ---------------------------------------

# limit to top N genes that give highest loadings...
N.genes <- 100

# get top N values from V. Therefore filter...
# Value calculated from Modulus of elements in vector to take into account
# both Re and Im

for (PCA in 1:10){
  top.N <- GetTopNValues(Mod(dat.pca$u[, PCA]), N.genes)  # is it right to get the Mod for TopNValues?
  
  # filter expression matrix for only top contributing genes
  # recreate approximation of gene exprs with PCA
  # Create this by outer product of eigengenes and eigensamples multiplied by eigenvalue
  eigenvalue <- dat.pca$d[PCA]
  eigengene <- dat.pca$v[, PCA]  # columns are eigengenes (genes are mix of these). v is 12x12
  # eigensample is large (k by 12), filter to top top.N genes...
  eigensample <- dat.pca$u[top.N$i, PCA]  # columns are eigensamples (samples are mix of these). u is top.N$i x 12
  # get pca.matrix by outer product of eigengene and eigensample. 
  pca.matrix <- eigenvalue * outer(eigensample, eigengene)
  colnames(pca.matrix) <- dat.tissuenames
  
  # order by phase angle
  order.phase <- order(Arg(eigensample))
  y <- 1:length(eigengene)  # length of 12, genes are mix of these.
  x <- 1:length(eigensample)  # length of top.N. samples are mix of these.
  image(x, y, Arg(pca.matrix[order.phase, ]), 
        main=paste('Phase angles: PCA:', PCA), 
        axes=FALSE, xlab="", ylab="")
  axis(1, at=x, labels=FALSE, tick=FALSE)
  axis(2, at=y, labels=FALSE, tick=FALSE)
  # Now draw the textual axis labels
  # gene labels
  text(x, par("usr")[3] - 1,
       labels=names(eigensample),
       srt=90, 
       pos=4,
       offset=0,
       xpd=TRUE, 
       cex=0.9) 
  # sample labels
  text(par("usr")[1] - 5, y, 
       labels = names(eigengene), 
       srt=0, 
       pos=3, 
       offset=0,
       xpd=TRUE, 
       cex=1.5) 
}

# # Clock Bmal
# pca1.arg <- t(dat.pca$u[, 1] %*% d.v)
# colnames(pca1.arg) <- dat.tissuenames
Bmal1 <- pca1.arg["Arntl", ]
Clock <- pca1.arg["Clock", ]

# Per1 Per2 Cry1
Per1 <-  pca1.arg["Per1", ]
Per2 <-  pca1.arg["Per2", ]
Cry1 <- pca1.arg["Cry1", ]

plot(c(Bmal1, Per1), col=c(rep(c('red', 'blue'), each=12)), xlim=c(-0.2, 0.2), ylim=c(-1, 1))
abline(v=0)
abline(h=0)
# 



