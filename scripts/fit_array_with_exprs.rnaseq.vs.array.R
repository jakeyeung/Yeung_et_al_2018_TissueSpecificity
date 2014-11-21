# Jake Yeung
# November 20 2014
# Linear fit


# Define constants --------------------------------------------------------

N.TISSUES <- 12
RNA.SEQ.INTERVAL <- 8
ARRAY.INTERVAL <- 2

# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'GetTissueNames.R'))
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names

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
scatter.outpath <- file.path(plot_dir, "scatter.rna.array.rnaseq.vs.array.pdf")
scatter.replicates.outpath <- file.path(plot_dir, "scatter.replicates.rnaseq.vs.array.pdf")
clock.genes.outpath <- file.path(plot_dir, "clock.genes.outpath.rnaseq.vs.array.pdf")
tissue.genes.outpath <- file.path(plot_dir, "tissue.genes.outpath.rnaseq.vs.array.pdf")


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

tissue.names <- sapply(colnames(array.exprs), function(x){
  tissue <- substring(x, 1, nchar(x) - 2)
})
tissue.names <- unique(unname(tissue.names))


# Get array times --------------------------------------------------------

array.times <- sapply(colnames(array.exprs), function(x){
  time <- substring(x, nchar(x) - 1, nchar(x))
  return(time)
})
(array.times <- unique(array.times))

# Get tissue names: RNASeq ------------------------------------------------

tissue.names.rnaseq <- sapply(colnames(rna.seq.exprs), function(x){
  tissue <- substring(x, 1, nchar(x) - 5)
  return(tissue)
})
tissue.names.rnaseq <- unique(tissue.names.rnaseq)


# Get RNASeq times --------------------------------------------------------

rna.seq.times <- unname(sapply(colnames(rna.seq.exprs), function(x){
  substr(x, nchar(x)-1, nchar(x))
}))
rna.seq.times <- unique(rna.seq.times)

# create greppable expression for times
rna.seq.times <- paste0('*', rna.seq.times, collapse='|')


# Subset microarray for rna.seq.times -------------------------------------

array.subset <- array.exprs[, grepl(rna.seq.times, colnames(array.exprs))]
Peek(array.subset)
Peek(rna.seq.exprs)

# Fit linear model --------------------------------------------------------

# get intersection of genes
common.genes <- intersect(rownames(array.subset), rownames(rna.seq.exprs))

array.subset.common.g <- as.matrix(array.subset[common.genes, ])
rna.seq.exprs.common.g <- as.matrix(rna.seq.exprs[common.genes, ])
array.common.g <- as.matrix(array.exprs[common.genes, ])  # for interpolation later.

coeff.mat <- matrix(NA, nrow=nrow(array.subset.common.g), ncol=2 * N.TISSUES)
dim(coeff.mat)

# set up row and colnames of coeff.mat
rownames(coeff.mat) <- rownames(array.subset.common.g)
coeff.mat.colnames <- rep(NA, 2 * N.TISSUES)
for (i in 1:N.TISSUES){
  tissue <- tissue.names[i]
  coeff.mat.colnames[i * 2] <- paste0(tissue, '_slope')
  coeff.mat.colnames[i * 2 - 1] <- paste0(tissue, '_intercept')
}
colnames(coeff.mat) <- coeff.mat.colnames

# rownames(array.subset.common.g)
for (gene.i in 1:nrow(array.subset.common.g)){
  # for each gene, do a fit for every tissue
  for (tissue.i in 1:N.TISSUES){
    tissue.i.end <- tissue.i * RNA.SEQ.INTERVAL
    tissue.i.start <- tissue.i.end - RNA.SEQ.INTERVAL + 1
    fit <- lm(rna.seq.exprs.common.g[gene.i, tissue.i.start:tissue.i.end] ~ array.subset.common.g[gene.i, tissue.i.start:tissue.i.end])
    intercept <- fit$coefficient[1]
    slope <- fit$coefficient[2]
    coeff.mat[gene.i, tissue.i * 2] <- slope
    coeff.mat[gene.i, tissue.i * 2 - 1] <- intercept
  }
}
Peek(coeff.mat)



# Adjust all microarray to RNAseq -----------------------------------------

array.exprs.adjusted <- matrix(NA, nrow=nrow(array.common.g), ncol=ncol(array.common.g),
                               dimnames=list(rownames(array.common.g),
                                             colnames(array.common.g)))
  
for (i in 1:N.TISSUES){
  tissue <- tissue.names[i]  # for coeff.mat and array
  # tissue.rnaseq <- tissue.names.rnaseq[i]  # for rnaseq
  intercept <- coeff.mat[, paste0(tissue, '_intercept')]
  slope <- coeff.mat[, paste0(tissue, '_slope')]
  tissue.exprs.array <- array.common.g[, grepl(paste0(tissue, '*'), 
                                               colnames(array.common.g))]
  tissue.exprs.array.normalized <- intercept + slope * tissue.exprs.array
  # write to adjusted exprs
  array.exprs.adjusted[, grepl(paste0(tissue, '*'), 
                                      colnames(array.exprs.adjusted))] <- tissue.exprs.array.normalized
#   tissue.exprs.rnaseq <- rna.seq.exprs.common.g[, grepl(paste0(tissue, '*'), 
#                                                         colnames(rna.seq.exprs.common.g))]
}
Peek(array.exprs.adjusted)


# Plot some genes for checking fit... -------------------------------------

# Look at negative guys!
negs <- which(array.exprs.adjusted < -18200, arr.ind = TRUE)
genes.i <- unique(negs[, "row"])
length(genes.i)

# plot for one tissue only
t <- 'Liver'
t.grep <- paste0(t, '*')

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
#   plot(array.subset.common.g[gene, grepl(t.grep, colnames(array.subset.common.g))], 
#        rna.seq.exprs.common.g[gene, grepl(t.grep, colnames(array.subset.common.g))], 
#        #        xlim=c(0, max(array.subset.common.g[gene, 1:24])),
#        #        ylim=c(0, max(rna.seq.exprs.common.g[gene, 1:24])),
#        main=paste(gene, 'Intercept=', intercept, 'Slope=', slope))
#   gfit <- lm(rna.seq.exprs.common.g[gene, grepl(t.grep, colnames(array.subset.common.g))] ~ array.subset.common.g[gene, grepl(t.grep, colnames(array.subset.common.g))])
#   abline(intercept, slope, lty='dotted')
#   abline(gfit)
  x <- array.exprs[gene, grepl(t.grep, colnames(array.exprs.subset))]
  print(length(x))
  y <- slope * x + intercept
  print(min(y))
  plot(x, y, main=paste("normal", gene, t))
  abline(h=0)
  
# log 2 transform that shit
  y[which(y < 1)] <- 1 
  plot(log2(x), log2(y), main=paste("log2", gene, t))
  # if you see two lines, one dotted and one solid, it means you did the fit wrong! They should be the same line.
}

# random genes for Adr
# set.seed(0)
random.genes <- sample(rownames(array.subset.common.g), 10)

# add clock genes
random.genes <- c(random.genes, 'Arntl', 'Per2', 'Per3', 'Dbp')

# plot for one tissue only
t <- 'Adr'
t.grep <- paste0(t, '*')
for (gene in random.genes){
  intercept <- coeff.mat[gene, paste0(t, '_intercept')]
  slope <- coeff.mat[gene, paste0(t, '_slope')]
  plot(array.subset.common.g[gene, grepl(t.grep, colnames(array.subset.common.g))], 
       rna.seq.exprs.common.g[gene, grepl(t.grep, colnames(array.subset.common.g))], 
       #        xlim=c(0, max(array.subset.common.g[gene, 1:24])),
       #        ylim=c(0, max(rna.seq.exprs.common.g[gene, 1:24])),
       main=paste(gene, 'Intercept=', intercept, 'Slope=', slope))
  gfit <- lm(rna.seq.exprs.common.g[gene, grepl(t.grep, colnames(array.subset.common.g))] ~ array.subset.common.g[gene, grepl(t.grep, colnames(array.subset.common.g))])
  abline(intercept, slope, lty='dotted')
  abline(gfit)
  # if you see two lines, one dotted and one solid, it means you did the fit wrong! They should be the same line.
}



# Scatterplot dat shit ----------------------------------------------------

genes.subset <- sample(rownames(array.exprs.adjusted), 0.1 * nrow(array.exprs.adjusted))

pdf('plots/adjust_scatter.pdf')
par(mfrow=c(1, 1))
# plot before
x1 <- unlist(array.exprs[genes.subset, grepl(rna.seq.times, colnames(array.exprs.adjusted))])
y <- unlist(rna.seq.exprs.common.g[genes.subset, grepl(rna.seq.times, colnames(rna.seq.exprs.common.g))])
plot(log2(x1), log2(y), xlab='Microarray unadjusted (log2)', ylab='RNA-Seq (log2)', main='Scatter plot: before adjustment')
x <- unlist(array.exprs.adjusted[genes.subset, grepl(rna.seq.times, colnames(array.exprs.adjusted))])
plot(log2(x), log2(y), xlab='Microarray adjusted (log2)', ylab='RNA-Seq (log2)', main='Scatter plot: after adjustment')
dev.off()


# Plot replicates ---------------------------------------------------------

# only loop til 50, sincei ts matching pair is 64
pdf('plots/adjust_exprs_pairs.pdf')
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
  # plot rnaseq
  t1.rnaseq <- 22
  t2.rnaseq <- 46
  x.rnaseq <- unlist(rna.seq.exprs.common.g[genes.subset, grepl(paste0('*', t1.rnaseq), colnames(rna.seq.exprs.common.g))])
  y.rnaseq <- unlist(rna.seq.exprs.common.g[genes.subset, grepl(paste0('*', t2.rnaseq), colnames(rna.seq.exprs.common.g))])
  plot(log2(x.rnaseq), log2(y.rnaseq), xlab=paste(t1.rnaseq, 'CT time'), ylab=paste(t2.rnaseq, 'CT time'), main='RNASeq (log2 DESeq-norm counts)')
}
dev.off()


# Plot clock genes --------------------------------------------------------

pdf('plots/before_after_adjusts_clock_genes.pdf')
par(mfrow=c(3,1))
clockgenes <- c('Dbp', 'Arntl', 'Clock', 'Per2', 'Per3', 'Cry2', 'Cry1')
for (gene in clockgenes){
  plot(log2(array.exprs[gene, ]), main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24))
  plot(log2(array.exprs.adjusted[gene, ]), main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24))
  plot(log2(rna.seq.exprs.common.g[gene, ]), main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8))
}
par(mfrow=c(1,1))
dev.off()

# Plot tissue-specific genes ----------------------------------------------

# Plot clock genes --------------------------------------------------------

pdf(clock.genes.outpath)
par(mfrow=c(3,1))
log2.array.exprs.adjusted <- array.exprs.adjusted
log2.array.exprs.adjusted[which(log2.array.exprs.adjusted < 0)] <- 0
log2.array.exprs.adjusted <- log2(log2.array.exprs.adjusted + 1)
clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12', 
                'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')
for (gene in clockgenes){
  plot(log2(array.exprs[gene, ]), main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, 14))
  #   plot(log2(array.exprs.adjusted[gene, ]), main=paste(gene, 'log2 exprs: array after adjustment'),
  #        col=rep(1:N.TISSUES, each=24), type='b')
  plot(log2.array.exprs.adjusted[gene, ], main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, 14))
  plot(log2(rna.seq.exprs.common.g[gene, ] + 1), main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8), type='b', ylim=c(0, 14))
}
par(mfrow=c(1,1))
dev.off()

# Plot tissue-specific genes ----------------------------------------------

pdf(tissue.genes.outpath)
par(mfrow=c(3,1))
tissue.genes <- c("Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
tissue.genes <- c(tissue.genes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R
for (gene in tissue.genes){
  plot(log2(array.exprs[gene, ]), main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24))
  plot(log2(array.exprs.adjusted[gene, ]), main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24))
  plot(log2(rna.seq.exprs.common.g[gene, ]), main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8))
}
dev.off()

