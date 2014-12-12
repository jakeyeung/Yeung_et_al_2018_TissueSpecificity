# Jake Yeung
# Dec 9 2014
# measure.entropy.loadings.R
# Measure entropy of left singular value modules to determine "tissue specificity" of modules.


# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'RemoveProblemGenes.R'))  # for getting tissue times for adjusted data
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'GetTissueTimes.R'))  # for getting tissue times for adjusted data
source(file.path(functions.dir, 'FourierFunctions.R'))  # for Fourier stuff
source(file.path(functions.dir, 'GetTopNVAlues.R'))  # get top genes in singular values


# Define constants --------------------------------------------------------

N.TIMEPTS <- 24  # 24 time points, 2 hours per time point (over 48 hrs)
N.TISSUES <- 12  # 12 tissues
PERIOD <- 24  # 24 hours
INTERVAL <- 2  # hours between samples
epsilon <- 1  # RNA-Seq log2 transform to prevent infinities

# Load data ---------------------------------------------------------------

# define dirs
data.dir <- "data"
array.fname <- "array.adj.0.07.txt"
rnaseq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"


# handle array data -------------------------------------------------------

array.dat <- read.table(file.path(data.dir, array.fname))
array.dat <- RemoveProblemGenes(array.dat)

# log2 transform
array.dat <- log2(array.dat + 1)

Peek(array.dat)
dat.tissuenames <- GetTissues(samp.names=colnames(array.dat))


# Handle RNA Seq data -----------------------------------------------------

rnaseq.dat <- read.table(file.path(data.dir, rnaseq.fname), header=TRUE, sep='\t')
rownames(rnaseq.dat) <- make.names(rnaseq.dat$gene, unique=TRUE)

drop.cols <- c("gene")
rnaseq.dat <- rnaseq.dat[, !(names(rnaseq.dat) %in% drop.cols)]

Peek(rnaseq.dat)  # expect gene names as row names, tissues in columns


# Project onto flat and rhythmic time components --------------------------

# RHYTHMIC
OMEGA <- 2 * pi / PERIOD
# OMEGA <- 0
dat.time.projected <- ProjectToPeriodicTime(as.matrix(array.dat), 
                                            N.TISSUES, 
                                            N.TIMEPTS, 
                                            INTERVAL, 
                                            OMEGA, 
                                            dat.tissuenames)


# SVD decomposition -------------------------------------------------------

dat.svd <- svd(Mod(dat.time.projected))  # not scaled
plot(dat.svd$d^2, type='o')  # Manual screeplot. 
print(lapply(dat.svd, dim))  # $u is 12x12, $v is k-genes x 12


# Get column and row names ------------------------------------------------

rownames(dat.svd$u) <- rownames(dat.time.projected)
colnames(dat.svd$v) <- colnames(dat.time.projected)


# Get top 100 genes -------------------------------------------------------

pdf("plots/entropy.gene.loadings.pdf")
top.n <- 100
components <- seq(1, length(dat.svd$d))
for (component in components){
  tops <- GetTopNValues(Mod(dat.svd$u[, component]), top.n)
  top.genes <- names(dat.svd$u[names(tops$vals), component])
  
  # Calculate entropy measure -----------------------------------------------
  # get "relative abundance" measurement for each gene across eigensamples
  p.total <- apply(dat.svd$u, 1, function(M){
    p.gene <- sum(Mod(M) ^ 2)
    return(p.gene)
  })
  
  # get matrix normalized by p.total
  p <- Mod(dat.svd$u) ^ 2 / p.total
  p <- p[complete.cases(p), ]
  
  H <- apply(p, 1, function(x){
    # Calculate entropy
    sum(x * log(1 / x)) 
  })
  
  par(mfrow = c(2, 1))
  plot(density(H), xlim=c(0, 2.5), main=paste("Density of H across all genes"))
  abline(v = mean(H))
  abline(v = median(H), col='blue')
  
  plot(density(H[top.genes]), xlim=c(0, 2.5), main=paste("Density of H across top", top.n, "genes. Component", component))
  abline(v = mean(H[top.genes]))
  abline(v = median(H[top.genes]), col='blue')
  
  dev.off() 
}