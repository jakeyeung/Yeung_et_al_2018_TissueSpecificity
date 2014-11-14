# Jake Yeung
# tissue_time_interactions_svd.R
# 2014-11-07
# Fix periodicity in time component. SVD on tissues.


# Constants ---------------------------------------------------------------

epsilon <- 1  # RNA-Seq log2 transform to prevent infinities

# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'PcaPlotFunctions.R'))  # for PCA and periodogram calcs
source(file.path(functions.dir, 'FourierFunctions.R'))  # for Fourier stuff
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
# PhaseToHsv package, loaded from github
# install_github("naef-lab/PhaseHSV")
library(PhaseHSV)

PlotComplex <- function(complex.matrix, gene.list, labels,
                        col="HSV", axis.min=-1.1, axis.max=1.1, 
                        main='Plot title'){
  # Plot genes on complex plane
  # 
  # ARGS:
  # complex matrix Gene expression. Genes are rows, samples are columns.
  #   expect rownames in complex matrix whcih are gene names (and match genelist)
  # gene.list Optionally filter by gene list
  # colors Colors, if "HSV" compute HSV from angles using PhaseToHsv
  
  if (missing(gene.list)){
    dat <- complex.matrix  
  } else {
    dat <- complex.matrix[gene.list, ]
  }
  if (missing(labels)){
    text.labels <- rownames(dat)  
  } else {
    text.labels <- labels
  }
  
  plot.colors <- hsv(h=PhaseToHsv(Arg(dat), -pi, pi), s=1, v=1)
    
  plot(dat, col=plot.colors, 
       xlim=c(axis.min, axis.max), 
       ylim=c(axis.min, axis.max), 
       pch=20,
       main=main)
  text(dat, 
       labels=text.labels, 
       pos=3)
  abline(v=0)
  abline(h=0)
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

PlotExprs <- function(cond.time.vec, 
                      N.TIMEPTS, 
                      INTERVAL,
                      N.CONDS, 
                      labels, 
                      main='Plot title'){
  # Plot expression across tissues over time.
  # 
  # ARGS:
  # cond.time.vec Vector of multiple conditions over time.
  # N.TIMEPTS Number of timepoints per tissue
  # N.CONDS Number of tissues
  # labels Tissue names
  # main Plot title
  # 
  # RETURN:
  # plot of exprs
  
  # create matrix from cond.time.vec
  exprs.mat <- matrix(cond.time.vec, 
                          ncol=N.TIMEPTS, 
                          nrow=N.CONDS, 
                          byrow=TRUE)
  time.mat <- matrix(seq(from=0, by=INTERVAL, length.out=N.TIMEPTS), nrow=N.TIMEPTS, ncol=N.CONDS)
  matplot(time.mat, t(exprs.mat), type=c("b"), pch=c(1:N.CONDS), col=1:N.CONDS, main=main,
          xlab='Time', ylab='Exprs')
  legend("topleft", legend=labels, col=1:N.CONDS, pch=c(1:N.CONDS)) # optional legend
}

GetTissueNames <- function(tissue.names.all){
  # Convert gene names Adr_CT22 to Adr for all tissue.names.all
  # 
  # Args:
  # tissue.names.all List of names like Adr_CT22 ...
  # 
  # Returns:
  # tissue.names List of names like Adr ...
  
  tissues.names <- strsplit(tissue.names.all, '_')
  tissues.names <- unique(unlist(lapply(tissues.names, '[[', 1)))
  return(tissues.names)
} 

# MAIN --------------------------------------------------------------------


# Define constants --------------------------------------------------------

N.TIMEPTS <- 24  # 24 time points, 2 hours per time point (over 48 hrs)
N.TISSUES <- 12  # 12 tissues
PERIOD <- 24  # 24 hours
INTERVAL <- 2  # hours between samples


# Load data ---------------------------------------------------------------

# define dirs
data_dir <- "data"
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
# OMEGA <- 0
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


# Plot MATRIX, TISSUE and GENE Loadings -----------------------------------

par(mfrow=c(1,1))
for (PCA in 1:5){
  # --------- BEGIN: GET EIGEN-MATRICES ----------------- # 
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

  # set colors, HUES, SATURATIOn, VALUE (HSV)
  # hues <- PhaseToHsv(Arg(pca.matrix), -pi, pi)
  hues <- seq(from=0, to=1, length.out=100)
  mycolors <- hsv(h=hues, s=1, v=1)
  
  y <- 1:length(eigengene)  # length of 12, genes are mix of these.
  x <- 1:length(eigensample)  # length of top.N. samples are mix of these.
  # --------- END: GET EIGEN-MATRICES ----------------- #
  
  # --------- BEGIN: PLOT MATRIX OF PHASE ANGLES ----------------- # 
  # 
  # order by phase angle
  order.phase <- order(Arg(eigensample))
  image(x, y, Arg(pca.matrix[order.phase, ]),
        col=mycolors,
        main=paste('Phase angles: PCA:', PCA), 
        axes=FALSE, xlab="", ylab="")
  axis(1, at=x, labels=FALSE, tick=FALSE)
  axis(2, at=y, labels=FALSE, tick=FALSE)
  # Now draw the textual axis labels
  # gene labels
  text(x, par("usr")[3] - 1.3,
       labels=names(eigensample)[order.phase],
       srt=90, 
       pos=4,
       offset=0,
       xpd=TRUE, 
       cex=0.9) 
  # sample labels
  text(par("usr")[1] - 3, y, 
       labels = names(eigengene), 
       srt=0, 
       pos=2, 
       offset=0,
       xpd=TRUE, 
       cex=1.5) 
  # ---------- END: PLOT MATRIX OF PHASE ANGLES ------------------- # 
  
  # ---- BEGIN: PLOT TISSUE AND GENE LOADINGS ON COMPLEX PLANE ---- # 
  # tissue loadings (eigengenes)
  PlotComplex(eigengene, axis.min=-0.5, axis.max=0.5, ,
              labels=dat.tissuenames,
              main=paste('Tissue loadings PCA:', PCA))
  
  # gene loadings (eigensamples)
  PlotComplex(eigensample, axis.min=-0.5, axis.max=0.5, 
              labels=names(eigensample),
              main=paste('Gene loadings PCA:', PCA))
  
  # ----- END: PLOT TISSUE AND GENE LOADINGS ON COMPLEX PLANE ----- # 
}


# Load RNA-Seq data -------------------------------------------------------

# Load data ---------------------------------------------------------------

# define dirs
fname.rnaseq <- "exprs_combined.genenames.txt"  # has duplicate gene names

# load data
data_path.rnaseq <- file.path(data_dir, fname.rnaseq)
print(paste("Reading data from,", data_path.rnaseq, "May take a few a minutes."))
dat.rnaseq <- read.table(data_path.rnaseq, header=TRUE)
print("Read data to memory.")



# If needed handle duplicate row names ------------------------------------

rownames(dat.rnaseq) <- make.names(dat.rnaseq$Gene.ID, unique=TRUE)

drop.cols <- c("Gene.ID")
dat.rnaseq <- dat.rnaseq[, !(names(dat.rnaseq) %in% drop.cols)]

Peek(dat.rnaseq)  # expect gene names as row names, tissues in columns


# log2 normalize ----------------------------------------------------------

dat.rnaseq <- log2(dat.rnaseq + epsilon)


# Plot clock gene exprs for microarray and RNA-Seq. -----------------------
# In RNASeq, BFat and BStm need to be swapped so it matches order with microarray
# easier for plotting visualization.

top.N <- GetTopNValues(Mod(dat.pca$u[, 1]), 20)
top.genes <- names(top.N$vals)

# filter original exprs data by top.N genes
dat.top.genes <- dat[top.genes, ]
dat.rnaseq.top.genes <- dat.rnaseq[top.genes, ]
  
# define tissues
tissues <- paste0(dat.tissuenames, '*')  # add * because we grep later
tissues.rnaseq <- paste0(GetTissueNames(colnames(dat.rnaseq)), '*')
# swap BFat with BStm for rnaseq. BFat at 3rd index. Bstm at 4th index.
tissues.rnaseq[c(3, 4)] <- tissues.rnaseq[c(4, 3)]
tissues.rnaseq

# grep tissues like Adr*|Hrt*| ... | BS*
dat.top.genes.tissuefilt <- dat.top.genes[, 
                                          (grepl(paste(tissues, collapse="|"), 
                                                  colnames(dat.top.genes)))]
dat.rnaseq.top.genes.tissuefilt <- dat.rnaseq.top.genes[, 
                                                        (grepl(paste(tissues.rnaseq, collapse="|"), 
                                                                colnames(dat.rnaseq.top.genes)))]

# replace Bfat with Bstm for rnaseq
bfat.grepper <- paste0("BFat", '*')
bstm.grepper <- paste0("Bstm", '*')
bfat.indices <- which(grepl(bfat.grepper, colnames(dat.rnaseq.top.genes)))
bstm.indices <- which(grepl(bstm.grepper, colnames(dat.rnaseq.top.genes)))
# swap data (but columns don't get swapped)
dat.rnaseq.top.genes.tissuefilt[, c(bfat.indices, bstm.indices)] <- dat.rnaseq.top.genes.tissuefilt[, c(bstm.indices, bfat.indices)]
# swap colnames
colnames(dat.rnaseq.top.genes.tissuefilt)
colnames(dat.rnaseq.top.genes.tissuefilt)[c(bfat.indices, bstm.indices)] <- colnames(dat.rnaseq.top.genes.tissuefilt)[c(bstm.indices, bfat.indices)]
colnames(dat.rnaseq.top.genes.tissuefilt)

# plot only genes with exprs in bith microarray and rnaseq
genes.to.plot <- intersect(rownames(dat.top.genes.tissuefilt), 
                           rownames(dat.rnaseq.top.genes.tissuefilt))

for (gene in genes.to.plot){
  # plot side by side.
  par(mfrow=c(1, 2))
  PlotExprs(dat.top.genes.tissuefilt[gene, ], 
            N.TIMEPTS, 
            INTERVAL,
            length(tissues), 
            labels=tissues,
            main=paste('Microarray. Gene:', gene))
  PlotExprs(dat.rnaseq.top.genes.tissuefilt[gene, ], 
            N.TIMEPTS=8, 
            INTERVAL=6,
            length(tissues.rnaseq), 
            labels=tissues.rnaseq,
            main=paste('RNASeq. Gene:', gene))
}




