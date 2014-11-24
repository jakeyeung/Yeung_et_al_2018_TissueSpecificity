# Jake Yeung
# November 20 2014
# Linear fit


# Define constants --------------------------------------------------------

N.TISSUES <- 12
N.SAMP.PER.TISSUE.RNASEQ <- 8
ARRAY.INTERVAL <- 2

# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'RegressionFunctions.R'))
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'PlotFunctions.R'))

# define directories ------------------------------------------------------

# define dirs
data_dir <- "data"
fname.rna.seq <- "exprs_combined.genenames.txt"
# fname.array <- "hogenesch_2014_rma.genenames.txt"
fname.array <- "hogenesch_2014_rma.genenames.colnameordered.txt"
# fname.array <- "all.hogenesch.seq.background.normalized.genenames.cut.filter.txt"
# fname.array <- "all.hogenesch.background.nonormalize.genenames.txt"  # normalize MANUALLY by setting NORMALIZE to TRUE
NORMALIZE <- FALSE

# Define output directories -----------------------------------------------

plot_dir <- "plots"
scatter.outpath <- file.path(plot_dir, "scatter.rna.array.rnaseq.vs.array2.pdf")
scatter.replicates.outpath <- file.path(plot_dir, "scatter.replicates.rnaseq.vs.array2.pdf")
clock.genes.outpath <- file.path(plot_dir, "clock.genes.outpath.rnaseq.vs.array2.pdf")
tissue.genes.outpath <- file.path(plot_dir, "tissue.genes.outpath.rnaseq.vs.array2.pdf")
fit.normal.log2.outpath <- file.path(plot_dir, "normal.log.plot.check.near.zero2.pdf")
tissue.genes.check.outpath <- file.path(plot_dir, "normal.log.plot.check.near.zero.tissue2.pdf")

# load data: RNASeq and microarray ----------------------------------------


# load data: rnaseq
rna.seq.path <- file.path(data_dir, fname.rna.seq)
print(paste("Reading data from,", rna.seq.path, "May take a few a minutes."))
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')
print("Read data to memory.")

# load data: microarray
array.path <- file.path(data_dir, fname.array)
print(paste("Reading data from,", array.path, "May take a few a minutes."))
array.exprs <- read.table(array.path, header=TRUE, sep='\t')  # has duplicate rownames
print("Read data to memory.")
Peek(array.exprs)

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$Gene.ID, unique=TRUE)

drop.cols <- c("Gene.ID")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns

# Handle duplicate rownames: MICROARRAY -------------------------------------

# first column contained gene names
rownames(array.exprs) <- make.names(array.exprs$gene, unique=TRUE)
# genes <- array.exprs$gene
# coords <- array.exprs$coordinates
# probeid <- array.exprs$ID_REF

drop.cols <- c("gene")
array.exprs <- array.exprs[, !(names(array.exprs) %in% drop.cols)]
Peek(array.exprs)

# un-Log2 Transform ----------------------------------------------------------

array.exprs <- 2 ^ (array.exprs)
# count exprs already in normal scale

# Handle array exprs tissue names -----------------------------------------

colnames(array.exprs) <- ShortenSampNames(colnames(array.exprs), show="tissue.time")
Peek(array.exprs)


# Get tissue names --------------------------------------------------------

tissue.names <- GetTissueNames(colnames(array.exprs), dat.type="array")


# Get array times --------------------------------------------------------

array.times <- GetTimes(colnames(array.exprs))


# Get tissue names: RNASeq ------------------------------------------------

tissue.names.rnaseq <- GetTissueNames(colnames(rna.seq.exprs), dat.type="rna.seq")



# Get RNASeq times --------------------------------------------------------

rna.seq.times <- GetTimes(colnames(rna.seq.exprs))


# Subset microarray for rna.seq.times -------------------------------------

# create greppable expression for times
rna.seq.times <- paste0('*', rna.seq.times, collapse='|')

array.subset <- array.exprs[, grepl(rna.seq.times, colnames(array.exprs))]
Peek(array.subset)
Peek(rna.seq.exprs)

# Fit linear model --------------------------------------------------------

# get intersection of genes
common.genes <- intersect(rownames(array.subset), rownames(rna.seq.exprs))

array.subset.common.g <- as.matrix(array.subset[common.genes, ])
rna.seq.exprs.common.g <- as.matrix(rna.seq.exprs[common.genes, ])
array.common.g <- as.matrix(array.exprs[common.genes, ])  # for interpolation later.

coeff.mat <- LmGeneTissue(array.subset.common.g, rna.seq.exprs.common.g, 
                           row.names=rownames(array.subset.common.g), 
                           tissue.names, n.samps=N.SAMP.PER.TISSUE.RNASEQ, 
                           coeff.mat=coeff.mat)


# Adjust all microarray to RNAseq -----------------------------------------

array.exprs.adjusted <- AdjustArrayToRnaSeq(array.common.g, coeff.mat, tissue.names)

# Create log2 transform ---------------------------------------------------

log2.array.exprs.adjusted <- array.exprs.adjusted
log2.array.exprs.adjusted[which(log2.array.exprs.adjusted < 0)] <- 0
log2.array.exprs.adjusted <- log2(log2.array.exprs.adjusted + 1)

# Scatterplot dat shit ----------------------------------------------------

genes.subset <- sample(rownames(array.exprs.adjusted), 0.1 * nrow(array.exprs.adjusted))

pdf(scatter.outpath)
par(mfrow=c(1, 1))
# plot before
x1 <- unlist(array.exprs[genes.subset, grepl(rna.seq.times, colnames(array.exprs.adjusted))])
y <- unlist(rna.seq.exprs.common.g[genes.subset, grepl(rna.seq.times, colnames(rna.seq.exprs.common.g))])
plot(log2(x1), log2(y), xlab='Microarray unadjusted (log2)', ylab='RNA-Seq (log2)', main='Scatter plot: before adjustment')
# plot after
x <- unlist(array.exprs.adjusted[genes.subset, grepl(rna.seq.times, colnames(array.exprs.adjusted))])
plot(log2(x), log2(y), xlab='Microarray adjusted (log2)', ylab='RNA-Seq (log2)', main='Scatter plot: after adjustment')
dev.off()


# Plot replicates ---------------------------------------------------------

# only loop til 50, sincei ts matching pair is 64
pdf(scatter.replicates.outpath)
par(mfrow=c(3, 1))
for (i in 1:(length(array.times) / 2)){
  t1 <- array.times[i]
  t2 <- array.times[i + 24 / ARRAY.INTERVAL]
  # plot before as comparison
  x1 <- unlist(array.exprs[genes.subset, grepl(paste0('*', t1), colnames(array.exprs.adjusted))])
  y1 <- unlist(array.exprs[genes.subset, grepl(paste0('*', t2), colnames(array.exprs.adjusted))])
  plot(log2(x1), log2(y1), xlab=paste(t1, 'CT time'), ylab=paste(t2, 'CT time'), main='Before adjustment (log2)')
  # plot after
  x <- unlist(array.exprs.adjusted[genes.subset, grepl(paste0('*', t1), colnames(array.exprs.adjusted))])
  y <- unlist(array.exprs.adjusted[genes.subset, grepl(paste0('*', t2), colnames(array.exprs.adjusted))])
  plot(log2(x), log2(y), xlab=paste(t1, 'CT time'), ylab=paste(t2, 'CT time'), main='After adjustment (log2)')
  # plot rnaseq at 22 and 46 for comparison
  t1.rnaseq <- 22
  t2.rnaseq <- 46
  x.rnaseq <- unlist(rna.seq.exprs.common.g[genes.subset, grepl(paste0('*', t1.rnaseq), colnames(rna.seq.exprs.common.g))])
  y.rnaseq <- unlist(rna.seq.exprs.common.g[genes.subset, grepl(paste0('*', t2.rnaseq), colnames(rna.seq.exprs.common.g))])
  plot(log2(x.rnaseq), log2(y.rnaseq), xlab=paste(t1.rnaseq, 'CT time'), ylab=paste(t2.rnaseq, 'CT time'), main='RNASeq (log2 DESeq-norm counts)')
}
dev.off()



# Plot clock genes --------------------------------------------------------

pdf(clock.genes.outpath)
par(mfrow=c(3,1))
clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12', 
                'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')
for (gene in clockgenes){
  plot(log2(array.exprs[gene, ]), main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, 14), ylab="log2 exprs", xlab=paste(tissue.names, collapse=" "))
  plot(log2.array.exprs.adjusted[gene, ], main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, 14), ylab="log2 exprs", xlab=paste(tissue.names, collapse=" "))
  plot(log2(rna.seq.exprs.common.g[gene, ] + 1), main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8), type='b', ylim=c(0, 14), ylab="log2 exprs", xlab=paste(tissue.names, collapse=" "))
}
par(mfrow=c(1,1))
dev.off()

# Plot some genes for checking fit... -------------------------------------

pdf(fit.normal.log2.outpath)
par(mfrow=c(2, 2))
# plot for one tissue only
t <- 'Liver'
t.rnaseq <- 'Liv'
t.grep <- t

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12', 
                'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')

# for (i in genes.i){
#   gene <- rownames(array.exprs.adjusted)[i]
for (gene in clockgenes){
  intercept <- coeff.mat[gene, paste0(t, '_intercept')]
  slope <- coeff.mat[gene, paste0(t, '_slope')]
  x <- array.exprs[gene, grepl(t.grep, colnames(array.exprs))]
  # get Liver names that have rnaseq and ones without... use for plotting them with different symbols
  x.rna.seq.i <- colnames(array.exprs.subset[gene, grepl(t.grep, colnames(array.exprs.subset))])  # subset only ones that have mRNA matched (to see range in microarray)
  # convert tissue names to indices
  x.rna.seq.i <- names(x) %in% x.rna.seq.i  # logical True/False
  # create plot symbols '*' if in x.rna.seq.i, '+' otherwise
  plot.symbols <- sapply(x.rna.seq.i, function(x){
    if (x == TRUE){
      symbol <- 1
    } else {
      symbol <- 8
    }
  })
  y <- slope * x + intercept
  y[which(y < 1)] <- 1 
  
  plot(x, y, main=paste("o=samp w/ RNASeq+array", gene), pch=c(plot.symbols), xlab="Observed microarray (normal)", ylab="Predicted expression (normal)")
  abline(h=0)
  
  # log 2 transform that shit
  
  plot(log2(x), log2(y), main=paste("log2", gene, t), pch=c(plot.symbols), xlab="Observed microarray (log2)", ylab="Predicted expression (log2)")
  # if you see two lines, one dotted and one solid, it means you did the fit wrong! They should be the same line.
  
  # Plot dat original
  PlotRnaMicroarrayFit(t, gene, coeff.mat, array.subset.common.g, rna.seq.exprs.common.g, t.rnaseq)
  
  # Plot Microarray all of them
  
  exprs.rnaseq <- rna.seq.exprs.common.g[gene, grepl(t.rnaseq, colnames(rna.seq.exprs.common.g))]
  plot(log2(x[x.rna.seq.i] + 1), log2(exprs.rnaseq + 1), 
       main=paste0("log2 exprs: raw data"), 
       xlab="Microarray exprs (log2)", 
       ylab="RNASeq (log2)")
  # plot(x, pch=plot.symbols, main='Normal exprs cross samples', ylab="Normal microarray exprs", xlab="Sample Index")
}
dev.off()

# Plot tissue-specific genes ----------------------------------------------

pdf(tissue.genes.outpath)
par(mfrow=c(3,1))
tissue.genes <- c("Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
tissue.genes <- c(tissue.genes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R
for (gene in tissue.genes){
  plot(log2(array.exprs[gene, ]), main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24), ylim=c(0, 15))
  plot(log2.array.exprs.adjusted[gene, ], main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24), ylim=c(0, 15))
  plot(log2(rna.seq.exprs.common.g[gene, ] + 1), main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8), ylim=c(0, 15))
}
dev.off()

pdf(tissue.genes.check.outpath)
par(mfrow=c(2,1))
t <- 'BFAT'
t.grep <- 'BFAT*'
for (gene in tissue.genes){
  intercept <- coeff.mat[gene, paste0(t, '_intercept')]
  slope <- coeff.mat[gene, paste0(t, '_slope')]
  x <- array.exprs[gene, grepl(t.grep, colnames(array.exprs))]
  # get Liver names that have rnaseq and ones without... use for plotting them with different symbols
  x.rna.seq.i <- colnames(array.exprs.subset[gene, grepl(t.grep, colnames(array.exprs.subset))])  # subset only ones that have mRNA matched (to see range in microarray)
  # convert tissue names to indices
  x.rna.seq.i <- names(x) %in% x.rna.seq.i  # logical True/False
  # create plot symbols '*' if in x.rna.seq.i, '+' otherwise
  plot.symbols <- sapply(x.rna.seq.i, function(x){
    if (x == TRUE){
      symbol <- 8
    } else {
      symbol <- 1
    }
  })
  y <- slope * x + intercept
  plot(x, y, main=paste("normal", gene, t), pch=c(plot.symbols), xlab="Observed microarray", ylab="Predicted expression")
  abline(h=0)
  
  # log 2 transform that shit
  y[which(y < 1)] <- 1 
  plot(log2(x), log2(y), main=paste("log2", gene, t), pch=c(plot.symbols), xlab="Observed microarray", ylab="Predicted expression")
  # if you see two lines, one dotted and one solid, it means you did the fit wrong! They should be the same line.
}
dev.off()
