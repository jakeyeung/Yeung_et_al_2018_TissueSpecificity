# noise_model.R
# Jake Yeung
# November 26 2014
# Plot variance versus mean for each tissue and each gene.

library(ggplot2)
library(plyr)

# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names

# define directories ------------------------------------------------------

# define dirs
data.dir <- "data"
fname.rna.seq <- "rna_seq_deseq_counts_colnames_fixed.txt"
fname.array <- "array_exprs_colnames_fixed2.txt"

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


# Calculate mean and variance for each gene per tissue --------------------

N <- nrow(rna.seq.exprs) * length(tissue.names)  # one measurement for each gene for all tissues.
mean.var <- list(mean=matrix(NA, nrow=nrow(rna.seq.exprs), ncol=length(tissue.names),
                             dimnames = list(rownames(rna.seq.exprs), 
                                             tissue.names)), 
                 var=matrix(NA, nrow=nrow(rna.seq.exprs), ncol=length(tissue.names),
                            dimnames = list(rownames(rna.seq.exprs), tissue.names)))

for (j in 1:length(tissue.names)){
  tissue <- tissue.names[j]
  gene.tissue.exprs <- rna.seq.exprs[, grepl(tissue, colnames(rna.seq.exprs))]
  # calculate mean and var, by row
  exprs.mean <- apply(gene.tissue.exprs, 1, mean)
  exprs.var <- apply(gene.tissue.exprs, 1, var)
  # append to matrix 
  mean.var$mean[, j] <- exprs.mean
  mean.var$var[, j] <- exprs.var
}
str(mean.var)

# Plot ggplot2 ------------------------------------------------------------

ggplot(mean.var.df, aes(x=mean, y=var)) + 
  geom_point(alpha=0.005) + 
  scale_x_log10() + 
  scale_y_log10()

# Make dataframe for ggplot2 ----------------------------------------------

mean.var.df <- data.frame(mean=as.vector(mean.var$mean), var=as.vector(mean.var$var))


# Bin and fit loess -------------------------------------------------------

n.per.bin <- 2000

mean.var.df <- mean.var.df[with(mean.var.df, order(mean)), ]
mean.var.df$bins.order <- factor(round_any(seq(1:nrow(mean.var.df)), 1000))
# bins.order <- factor(floor(order(mean.var.df$mean) / n.per.bin))
# mean.var.df <- cbind(mean.var.df, bins.order)

# Bin mean and bin variance: by order -------------------------------------

prob='mean'
bin.mean <- with(mean.var.df, (tapply(mean, bins.order, function(x){
   mean.quantiles <- median(x)
})))
bin.var <- with(mean.var.df, (tapply(var, bins.order, function(x){
  var.quantiles <- median(x)
  return(var.quantiles)
})))

# remove last point
bin.mean <- bin.mean[1:(length(bin.mean) - 1)]
bin.var <- bin.var[1:(length(bin.var) - 1)]

plot(bin.mean, bin.var, 
     main=paste0(n.per.bin, " genes per bin. ", prob, " quantiles."))

m0 <- 1
b0 <- 0.1 
k0 <- 0.1
m.min <- 0
b.min <- -Inf
k.min <- 0
m.max <- Inf
b.max <- Inf
k.max <- Inf
# Fit with constraint
f.parab <- function(x, bin.mean) x[1] * ( bin.mean - x[2] ) ^ 2 + x[3]
S.parab <- function(x) sum((bin.var - f(x, bin.mean)) ^ 2)
fit <- optim(c(m0, b0, k0), S.parab, method="L-BFGS-B", lower=c(m.min, b.min, k.min), upper=c(m.max, b.max, k.max))
str(fit)
x <- bin.mean
y <- f.parab(fit$par, bin.mean)
lines(x, y)     

# Get mean exprs across samples -------------------------------------------

# Create subset and common genes so each microarray point has corresponding rnaseq
common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))
common.samples <- intersect(colnames(array.exprs), colnames(rna.seq.exprs))
array.exprs.subset.common.g <- 2^(array.exprs[common.genes, common.samples])  # normal scale
rna.seq.exprs.common.g <- rna.seq.exprs[common.genes, ]
Peek(array.exprs.subset.common.g)
Peek(rna.seq.exprs.common.g)

# Take avg across samples
array.avg.exprs <- apply(array.exprs.subset.common.g, 1, mean)
rna.avg.exprs <- apply(rna.seq.exprs.common.g, 1, mean)


# Plot array and rna avg --------------------------------------------------

plot(array.avg.exprs, 
     rna.avg.exprs,
     pch=46, cex=2, log="xy")

# Fit noise model ---------------------------------------------------------

# model: R = a ( m - b ), b < min(m) in linear scale

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
                # 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
# clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')

R.hat <- function(m, a, b) a * (m - b)
for (gene in clockgenes){
  R <- rna.seq.exprs.common.g[gene, ]
  M <- array.exprs.subset.common.g[gene, ]
  M.full <- 2^array.exprs[gene, ]
  
  a.init <- 1.01  # expect to be > 1
  b.init <- 0.5 * min(M.full)  # background some fraction of min exprs
  
  a.min <- 0
  a.max <- 2^10
  b.min <- 0  # background can't be negative
  b.max <- min(M.full)  # background can't be larger than min exprs
  
  # calculate variance fom loess, giving it R. If outside of
  # interpolation range, set Ri = Rmin or Ri = Rmax
  # min.R <- 321
  # min.R <- min(fit$x)  # results in negative
  # max.R <- max(fit$x)
  
  # R.adj <- R
  # R.adj[which(R < min.R)] <- min.R
  # R.adj[which(R > max.R)] <- max.R
  
  # R.var <- predict(fit, data.frame(bin.mean=unlist(R.adj)))
  # R.var <- b0 + b1 * R + b2 * R ^ 2
  R.var <- f.parab(fit$par, R)

  S <- function(x) sum((R - R.hat(M, x[1], x[2])) ^ 2 / R.var)
  
  a.b <- optim(c(a.init, b.init), S, method="L-BFGS-B",
               lower=c(a.min, b.min),
               upper=c(a.max, b.max))
  a.hat <- a.b$par[1]
  b.hat <- a.b$par[2]
  conv <- a.b$convergence
  
  R.predict <- R.hat(M.full, a.hat, b.hat)
  
  plot(unlist(M), unlist(R), main=paste0("gene=", gene, 
                                         " a=", signif(a.hat, 2), 
                                         " b=", signif(b.hat, 2),
                                         " min.predict=", signif(min(R.predict), 2),
                                         " converge=", conv))
  lines(M.full, R.predict, lwd=3, col='red')
}
