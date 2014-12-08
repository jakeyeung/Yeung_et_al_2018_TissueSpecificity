# Jake Yeung
# linear fit: in log2 scale.
# Uses same functions as fit_with_exprs.rna.seq.vs.array.R

# Define constants --------------------------------------------------------

N.TISSUES <- 12
N.SAMP.PER.TISSUE.RNASEQ <- 8
ARRAY.INTERVAL <- 2
clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12', 
                'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')

tissue.genes <- c("Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
tissue.genes <- c(tissue.genes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R

# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'PlotFunctions.R'))

# define directories ------------------------------------------------------

# define dirs
data_dir <- "data"
fname.rna.seq <- "exprs_combined.genenames.txt"
fname.array <- "hogenesch_2014_rma.genenames.colnameordered.txt"

# Define output directories -----------------------------------------------

plot_dir <- "plots"
plot.tissue.array <- "Lung"
plot.tissue.rna.seq <- "Lun"
scatter.outpath <- file.path(plot_dir, "scatter.rna.array.rnaseq.vs.array.log2.fit.2.pdf")
scatter.replicates.outpath <- file.path(plot_dir, "scatter.replicates.rnaseq.vs.array2.log2.fit.2.pdf")
clock.genes.outpath <- file.path(plot_dir, "clock.genes.outpath.rnaseq.vs.array2.log2.fit.2.pdf")
tissue.genes.outpath <- file.path(plot_dir, "tissue.genes.outpath.rnaseq.vs.array2.log2.fit.2.pdf")
fit.normal.log2.outpath <- file.path(plot_dir, "normal.log.plot.check.near.zero2.log2.fit.lung.pdf")
tissue.genes.check.outpath <- file.path(plot_dir, "normal.log.plot.check.near.zero.tissue2.log2.fit.2.pdf")
rna.seq.vs.after.adj.overlay.outpath <- file.path(plot_dir, "rna.seq.vs.array.adj.log2.fit.2.pdf")
before.vs.after.adj.overlay.outpath <- file.path(plot_dir, "before.vs.array.adj.log2.fit.2.pdf")

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

#TODO: these duplicate gene names may be useful...

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


# Put into log2 scale -----------------------------------------------------

# microarray already log2 scale

rna.seq.exprs <- log2(rna.seq.exprs + 1)

# Fit linear model --------------------------------------------------------

# get intersection of genes
common.genes <- intersect(rownames(array.subset), rownames(rna.seq.exprs))

array.subset.common.g <- as.matrix(array.subset[common.genes, ])
rna.seq.exprs.common.g <- as.matrix(rna.seq.exprs[common.genes, ])
array.common.g <- as.matrix(array.exprs[common.genes, ])  # for interpolation later.

coeff.mat <- LmGeneTissue(array.subset.common.g, rna.seq.exprs.common.g, 
                          row.names=rownames(array.subset.common.g), 
                          tissue.names, n.samps=N.SAMP.PER.TISSUE.RNASEQ)

# Adjust all microarray to RNAseq -----------------------------------------

array.exprs.adjusted <- AdjustArrayToRnaSeq(array.common.g, coeff.mat, tissue.names)


# Plot clock genes --------------------------------------------------------

pdf(clock.genes.outpath)
par(mfrow=c(3,1))
for (gene in clockgenes){
  # Array before adjustment
  plot(t(as.matrix(array.exprs[gene, ])), main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, 14), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # Array after adjustment
  plot(as.matrix(array.exprs.adjusted[gene, ]), main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, 14), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # RNA Seq
  plot(as.matrix(rna.seq.exprs.common.g[gene, ]), main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8), type='b', ylim=c(0, 14), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
}
par(mfrow=c(1,1))
dev.off()


# Plot clock genes: overlap -----------------------------------------------

# plot rna-seq vs after adjustment

pdf(rna.seq.vs.after.adj.overlay.outpath)
for (gene in clockgenes){
  # RNA Seq
  rna.seq <- matrix(NA, nrow=1, ncol=ncol(array.exprs.adjusted), dimnames=list(gene, colnames(array.exprs.adjusted)))
  rna.seq[gene, colnames(array.exprs.subset.common.g)] <- as.matrix(rna.seq.exprs.common.g[gene, ])
  
  plot(seq(1:length(rna.seq)), rna.seq, main=paste(gene, 'black=rnaseq, red=array after adjust'),
       col=1, type='o', ylim=c(0, 14), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  lines(as.matrix(array.exprs.adjusted[gene, ]), col=2, pch=22, type='o')
}
dev.off()

# plot before vs after adjustment
pdf(before.vs.after.adj.overlay.outpath)
for (gene in clockgenes){
  # Array before adjustment
  plot(t(as.matrix(array.exprs[gene, ])), main=paste(gene, 'black=array before adjust, red=array after adjust'),
       col=1, type='b', ylim=c(0, 14), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  lines(as.matrix(array.exprs.adjusted[gene, ]), col=2, pch=22, type='o')
}
dev.off()


# PLOT DIAGNOSTIC PLOTS: clock genes ---------------------------------------------------

PlotFitDiagnostics(array.common.g, 
                   array.subset, 
                   rna.seq.exprs.common.g, 
                   clockgenes, 
                   coeff.mat,
                   outpath=fit.normal.log2.outpath, 
                   tissue=plot.tissue.array, 
                   tissue.rna.seq=plot.tissue.rna.seq,
                   allow.negs=TRUE)


# Plot diagnostic plots: tisues -------------------------------------------


pdf(tissue.genes.outpath)
par(mfrow=c(3,1))
for (gene in tissue.genes){
  plot(as.matrix(array.common.g[gene, ]), main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24), ylim=c(0, 15))
  plot(as.matrix(array.exprs.adjusted[gene, ]), main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24), ylim=c(0, 15))
  plot(t(as.matrix(rna.seq.exprs[gene, ])), main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8), ylim=c(0, 15))
}
dev.off()

PlotFitDiagnostics(array.common.g, 
                   array.subset, 
                   rna.seq.exprs.common.g, 
                   tissue.genes, 
                   coeff.mat,
                   outpath=tissue.genes.check.outpath, 
                   tissue="BFAT", 
                   tissue.rna.seq="BFat",
                   allow.negs=TRUE)
