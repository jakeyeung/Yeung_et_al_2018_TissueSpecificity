# noise_model.R
# Jake Yeung
# November 26 2014
# Plot variance versus mean for each tissue and each gene.

Rprof("profile1.out")

library(ggplot2)
library(plyr)
library(parallel)

# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'RegressionFunctions.R'))  # optim function

# define directories ------------------------------------------------------

# define dirs
data.dir <- "data"
fname.rna.seq <- "rna_seq_deseq_counts_colnames_fixed.txt"
fname.array <- "array_exprs_colnames_fixed.txt"
clock.genes.fit.outpath <- file.path("plots", "clock.genes.noise.model.fit.pdf")
clock.genes.outpath <- file.path("plots", "clock.genes.noise.model.before.after.pdf")
array.before.after.outpath <- file.path("plots", "array.before.after.noise.model.pdf")
rna.seq.array.outpath <- file.path("plots", "rna.seq.array.noise.model.pdf")
scatterplot.mean.variance.outpath <- file.path("plots", "scatterplot.mean.variance.noise.model.log2.pdf")
mean.var.fit.outpath <- file.path("plots", "scatterplot.mean.variance.parabola.fit.pdf")

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

# # remove last point
bin.mean <- bin.mean[1:(length(bin.mean) - 1)]
bin.var <- bin.var[1:(length(bin.var) - 1)]

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
S.parab <- function(x) sum((bin.var - f.parab(x, bin.mean)) ^ 2)
fit <- optim(c(m0, b0, k0), S.parab, method="L-BFGS-B", lower=c(m.min, b.min, k.min), upper=c(m.max, b.max, k.max))
str(fit)
x <- bin.mean
y <- f.parab(fit$par, bin.mean)

# plot dat
pdf(mean.var.fit.outpath)
plot(bin.mean, bin.var, 
     main=paste0('y=m(x-b)^2 + k. m=', signif(fit$par[1], 2), 
                 ' b=', signif(fit$par[2], 2), 
                 ' k=', signif(fit$par[3], 2)),
     xlim=c(0, 5000), ylim=c(0, 250000))
lines(x, y)     
dev.off()

# Get mean exprs across samples -------------------------------------------

# Create subset and common genes so each microarray point has corresponding rnaseq
common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))
common.samples <- intersect(colnames(array.exprs), colnames(rna.seq.exprs))
array.exprs.subset.common.g <- 2^(array.exprs[common.genes, common.samples])  # normal scale
rna.seq.exprs.common.g <- rna.seq.exprs[common.genes, ]
Peek(array.exprs.subset.common.g)
Peek(rna.seq.exprs.common.g)

# put all array exprs in normal scale
array.exprs <- 2 ^ array.exprs  # normal scale

# Take avg across samples
array.avg.exprs <- apply(array.exprs.subset.common.g, 1, mean)
rna.avg.exprs <- apply(rna.seq.exprs.common.g, 1, mean)


# Plot ggplot2 ------------------------------------------------------------

pdf(scatterplot.mean.variance.outpath)
ggplot(mean.var.df, aes(x=mean, y=var)) + 
  geom_point(alpha=0.005) + 
  scale_x_log10() + 
  scale_y_log10()
dev.off()


# Plot array and rna avg --------------------------------------------------

plot(array.avg.exprs, 
     rna.avg.exprs,
     pch=46, cex=2, log="xy")


# Clock genes --------------------------------------------------------

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')

# Fit noise model ---------------------------------------------------------

# model: R = a ( m - b ), b < min(m) in linear scale

# init coefficient matrix
coeff.mat <- data.frame(a.hat=rep(NA, length(common.genes)), 
                        b.hat=rep(NA, length(common.genes)), 
                        convergence=rep(NA, length(common.genes)), 
                        row.names=common.genes)
# init adjusted array
array.adj <- matrix(NA, nrow=length(common.genes), ncol=ncol(array.exprs), 
                    dimnames=list(common.genes, colnames(array.exprs)))

# array.adj <- mclapply(list(common.genes), ConstrainedFitWithNoise, 
#                            array.subset=array.exprs.subset.common.g, 
#                            array.full=array.exprs,
#                            rna.seq=rna.seq.exprs.common.g, 
#                            noise.function=f.parab, 
#                            noise.params=fit$par)
# array.adj <- array.adj[[1]]

R.hat <- function(m, a, b) a * (m - b) # RNA to Microarray model.
S <- function(x) sum((R - R.hat(M, x[1], x[2])) ^ 2 / R.var)  # optimization equation

count <- 1
for (gene in common.genes){
  R <- rna.seq.exprs.common.g[gene, ]
  M <- array.exprs.subset.common.g[gene, ]
  M.full <- array.exprs[gene, ]
  
  a.init <- 1.01  # expect to be > 1
  b.init <- 0.5 * min(M.full)  # background some fraction of min exprs
  
  a.min <- 0
  a.max <- 2^10
  b.min <- 0  # background can't be negative
  b.max <- min(M.full)  # background can't be larger than min exprs
  
  R.var <- f.parab(fit$par, R)
  
  a.b <- optim(c(a.init, b.init), S, method="L-BFGS-B",
               lower=c(a.min, b.min),
               upper=c(a.max, b.max))
  a.hat <- a.b$par[1]
  b.hat <- a.b$par[2]
  conv <- a.b$convergence
  
  R.predict <- R.hat(M.full, a.hat, b.hat)
  
  # add to coeff.mat
  coeff.mat[gene, "a.hat"] <- a.hat
  coeff.mat[gene, "b.hat"]  <- b.hat
  coeff.mat[gene, "convergence"] <- conv
  
  array.adj[gene, ] <- R.predict
  count <- count + 1
  if (count %% 1000 == 0){
    print(count)
  }
}

# save coeff.mat for later
save(coeff.mat, file="coeff.mat.noise.model.RData")


# Plot clock genes fit ----------------------------------------------------


pdf(clock.genes.fit.outpath)
for (gene in clockgenes){
  M <- array.exprs.subset.common.g[gene, ]
  R <- rna.seq.exprs.common.g[gene, ]
  a.hat <- coeff.mat[gene, "a.hat"]
  b.hat <- coeff.mat[gene, "b.hat"]
  M.full <- array.exprs[gene, ]
  R.predict <- R.hat(M.full, a.hat, b.hat)
  plot(unlist(M), unlist(R), main=paste0("gene=", gene, 
                                         " a=", signif(a.hat, 2), 
                                         " b=", signif(b.hat, 2),
                                         " min.predict=", signif(min(R.predict), 2)))
  lines(M.full, R.predict, lwd=3, col='red')
}
dev.off()


# Plot before after adjust and RNASeq: Log2 -------------------------------

N.TISSUES <- 12

pdf(clock.genes.outpath)
par(mfrow=c(3,1))
for (gene in clockgenes){
  # Array before adjustment
  plot(unlist(log2(array.exprs[gene, ])), main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, 14), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # Array after adjustment
  plot(unlist(log2(array.adj[gene, ])), main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, 14), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # RNA Seq
  plot(unlist(log2(rna.seq.exprs.common.g[gene, ] + 1)), main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8), type='b', ylim=c(0, 14), ylab="log2 exprs",
       xlab=paste(tissue.names, collapse=" "))
}
par(mfrow=c(1,1))
dev.off()


# Plot before and after on same scale -------------------------------------

pdf(array.before.after.outpath)
for (gene in clockgenes){
  # Array before adjustment
  plot(unlist(log2(array.exprs[gene, ])), main=paste(gene, 'black=array before, red=array after'),
       col=1, type='b', ylim=c(0, 14), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # Array after adjustment
  lines(unlist(log2(array.adj[gene, ])),
       col=2, pch=23, type='o')
}
dev.off()

pdf(rna.seq.array.outpath)
for (gene in clockgenes){
  # RNA Seq
  plot(unlist(log2(rna.seq.exprs.common.g[gene, ] + 1)), main=paste(gene, 'black=rnaseq, red=array after'),
       col=1, type='b', ylim=c(0, 14), ylab="log2 exprs",
       xlab=paste(tissue.names, collapse=" "))
  lines(unlist(log2(array.adj[gene, common.samples])),
        col=2, pch=23, type='o')
}
dev.off()

# pdf(rna.seq.array.outpath2)
# for (gene in clockgenes){
#   # RNA Seq
#   rna.seq <- matrix(NA, nrow=1, ncol=ncol(array.adj), dimnames=list(gene, colnames(array.adj)))
#   rna.seq[gene, colnames(array.adj)] <- as.matrix(rna.seq.exprs.common.g[gene, ])
#   
#   plot(seq(1:length(rna.seq)), unlist(log2(rna.seq + 1)), main=paste(gene, 'black=rnaseq, red=array after'),
#        col=1, type='b', ylim=c(0, 14), ylab="log2 exprs", 
#        xlab=paste(tissue.names, collapse=" "))
#   lines(unlist(log2(array.adj[gene, ])),
#         col=2, pch=23, type='o')
# }
# dev.off()



Rprof(NULL)
summaryRprof("profile1.out", lines = "show")



