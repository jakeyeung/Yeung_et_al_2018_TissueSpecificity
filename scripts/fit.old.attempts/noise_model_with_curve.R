# noise_model_with_curve.R
# Jake Yeung
# November 28 2014
# Plot variance versus mean for each tissue and each gene.
# include saturation in model of parameter...

# Clock genes --------------------------------------------------------

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')


# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'RegressionFunctions.R'))  # optim function
source(file.path(functions.dir, 'GeneTissueCalculations.R'))

# define directories ------------------------------------------------------

# define dirs
data.dir <- "data"
fname.rna.seq <- "rna_seq_deseq_counts_colnames_fixed.txt"
fname.array <- "array_exprs_colnames_fixed.txt"
clock.genes.fit.outpath <- file.path("plots", "clock.genes.noise.model.fit.curve.pdf")
clock.genes.outpath <- file.path("plots", "clock.genes.noise.model.before.after.curve.pdf")
array.before.after.outpath <- file.path("plots", "array.before.after.noise.model.curve.pdf")
rna.seq.array.outpath <- file.path("plots", "rna.seq.array.noise.model.curve.pdf")
scatterplot.mean.variance.outpath <- file.path("plots", "scatterplot.mean.variance.noise.model.log2.curve.pdf")
mean.var.fit.outpath <- file.path("plots", "scatterplot.mean.variance.parabola.fit.curve.pdf")

# Load file ---------------------------------------------------------------

# load data: rnaseq
rna.seq.path <- file.path(data.dir, fname.rna.seq)
print(paste("Reading data from,", rna.seq.path, "May take a few a minutes."))
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')
print("Read data to memory.")

array.path <- file.path(data.dir, fname.array)
array.exprs <- read.table(array.path, header=TRUE, sep='\t')

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns

# Handle duplicate rownames: MICROARRAY -------------------------------------

# first column contained gene names
rownames(array.exprs) <- make.names(array.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
array.exprs <- array.exprs[, !(names(array.exprs) %in% drop.cols)]
Peek(array.exprs)


# Get column names --------------------------------------------------------

tissue.names <- GetTissueNames(colnames(rna.seq.exprs)) 



# Get subsets of exprs ----------------------------------------------------

# Create subset and common genes so each microarray point has corresponding rnaseq
common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))
common.samples <- intersect(colnames(array.exprs), colnames(rna.seq.exprs))
array.exprs.subset.common.g <- 2^(array.exprs[common.genes, common.samples])  # normal scale
rna.seq.exprs.common.g <- rna.seq.exprs[common.genes, ]
Peek(array.exprs.subset.common.g)
Peek(rna.seq.exprs.common.g)

# Calculate MeanVariance by tissue ----------------------------------------

mean.var.df <- GetMeanVarByTissues(exprs=rna.seq.exprs, tissue.names)

# Bin and fit loess -------------------------------------------------------

n.per.bin <- 150

mean.var.df <- mean.var.df[with(mean.var.df, order(mean)), ]
mean.var.df$bins.order <- factor(round_any(seq(1:nrow(mean.var.df)), n.per.bin))

# Bin mean and bin variance: by order -------------------------------------

bin.mean <- with(mean.var.df, (tapply(mean, bins.order, function(x){
  mean.quantiles <- median(x)
})))
bin.var <- with(mean.var.df, (tapply(var, bins.order, function(x){
  var.quantiles <- median(x)
  return(var.quantiles)
})))


# Fit loess noise model ---------------------------------------------------

fit.noise <- loess(bin.var ~ bin.mean,
                   control=loess.control(surface="direct"))

# Plot fit ----------------------------------------------------------------


pdf(mean.var.fit.outpath)
x <- seq(min(rna.seq.exprs.common.g), max(rna.seq.exprs.common.g), 100)  # plot full range
y <- predict(fit.noise, x)
plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
     xlab="bin.mean", ylab="bin.var",
     xlim=c(0, 3000), ylim=c(0, 50000))
points(bin.mean, bin.var)
plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
     xlab="bin.mean", ylab="bin.var")
points(bin.mean, bin.var)
plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
     xlab="bin.mean", ylab="bin.var", log="xy")
points(bin.mean, bin.var)

dev.off()


# Fit curve ---------------------------------------------------------------

pdf(clock.genes.fit.outpath)
for (gene in clockgenes){
  b0 <- 16
  a0 <- 10000
  k0 <- 10000
  init.vals <- c(a0, b0, k0)
  # microarray rnaseq model. 
  m <- function(r, a, b, k) b + (a * r / k) / ( 1 + r / k)
  # define training set
  M <- array.exprs.subset.common.g[gene, ]
  R <- rna.seq.exprs.common.g[gene, ]
  # calcualte noise
  R.var <- predict(fit.noise, unlist(R))
  S.m <- function(x) sum(((M - m(R, x[1], x[2], x[3])) ^ 2) / R.var)
  fit.S <- optim(init.vals, S.m, method='L-BFGS-B')
  print(paste(gene, fit.S$convergence))
  plot(unlist(R), unlist(M), main=paste0(gene, 
                                         " m=b+(ar/k)/(1+r/k)", 
                                         " a=", signif(fit.S$par[1], 2), 
                                         " b=", signif(fit.S$par[2], 2), 
                                         " k=", signif(fit.S$par[3]), 2), 
       xlab="RNA-Seq", ylab="Microarray")
  lines(sort(unlist(R)), m(sort(unlist(R)), fit.S$par[1], fit.S$par[2], fit.S$par[3]))
}
dev.off()
