# Jake Yeung
# Nov 28 2014
# Fit linear model

library(plyr)
library(ggplot2)

# Clock genes --------------------------------------------------------

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')


# Tissue genes ------------------------------------------------------------

tissuegenes <- c("Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
tissuegenes <- c(tissue.genes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R


# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'RegressionFunctions.R'))  # optim function
source(file.path(functions.dir, 'GeneTissueCalculations.R'))
source(file.path(functions.dir, 'PlotFunctions.R'))

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



# log2 transform RNA-Seq --------------------------------------------------

rna.seq.exprs <- log2(rna.seq.exprs + 1)
# Create subset and common genes so each microarray point has corresponding rnaseq
common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))
common.samples <- intersect(colnames(array.exprs), colnames(rna.seq.exprs))
array.exprs.subset.common.g <- array.exprs[common.genes, common.samples]  # log2 scale
rna.seq.exprs.common.g <- rna.seq.exprs[common.genes, ]
Peek(array.exprs.subset.common.g)
Peek(rna.seq.exprs.common.g)


# Noise model: log2 scale -------------------------------------------------

mean.var.df <- GetMeanVarByTissues(exprs=rna.seq.exprs, tissue.names)


# Plot scatters and density -----------------------------------------------

# ggplot(mean.var.df, aes(x=mean)) + geom_density()
# ggplot(mean.var.df, aes(x=mean, y=var)) + geom_point(alpha=0.005)

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



# Linear regression -------------------------------------------------------

coeff.mat <- matrix(0, nrow=length(common.genes), ncol=2,
                    dimnames=list(common.genes, 
                                  c("slope", "intercept")))
for (gene in common.genes){
  R <- unlist(rna.seq.exprs.common.g[gene, ])
  M <- unlist(array.exprs.subset.common.g[gene, ])
  fit <- lm(R ~ M)
  coeff.mat[gene, "intercept"] <- fit$coefficients[1]
  coeff.mat[gene, "slope"] <- fit$coefficients[2]
}

coeff.mat[clockgenes, ]

# Visualize clock genes ---------------------------------------------------

# log(r) = a log(m) + b

for (gene in clockgenes){
  # create R vs M full 288, R = 0 for "missing" values... for plotting
  R.full <- matrix(0, nrow=1, ncol=ncol(array.exprs), 
              dimnames=list(gene, colnames(array.exprs)))
  R.full[gene, common.samples] <- as.matrix(rna.seq.exprs[gene, common.samples])
  M.full <- as.matrix(array.exprs[gene, ])
  # create M and R for fitting...
  R <- as.matrix(rna.seq.exprs[gene, common.samples])
  M <- as.matrix(array.exprs.subset.common.g[gene, ])
  # create plot symbols
  plot.symbols <- matrix(8, nrow=1, ncol=ncol(array.exprs),
                         dimnames=list(gene, colnames(array.exprs)))
  plot.symbols[gene, common.samples] <- 1
  plot(M.full, R.full, main=paste(gene, "Microarray vs RNA-Seq"), 
       xlab="Microarray (log2)", 
       ylab="RNA-Seq (log2)",
       pch=plot.symbols)
  int <- coeff.mat[gene, "intercept"]
  slope <- coeff.mat[gene, "slope"]
  abline(int, slope)
}


# Adjust my microarray ----------------------------------------------------

array.adj <- matrix(NA, nrow=length(common.genes), ncol=ncol(array.exprs), 
                    dimnames=list(common.genes, colnames(array.exprs)))
for (gene in common.genes){
  int <- coeff.mat[gene, "intercept"]
  slope <- coeff.mat[gene, "slope"]
  M.orig <- as.matrix(array.exprs[gene, ])
  M.adj <- slope * M.orig + int
  array.adj[gene, ] <- M.adj
}
save(array.adj, file='array.adj.lm.fit.log2.by.probe.RData')

# Plot clock genes --------------------------------------------------------

N.TISSUES <- 12
pdf('plots/clock.genes.log2.lin.fit.by.probe.pdf')
par(mfrow=c(3,1))
for (gene in clockgenes){
  PlotBeforeAfter(array.exprs, array.adj, rna.seq.exprs)
}
dev.off()

pdf('plots/tissue.genes.log2.lin.fit.by.probe.pdf')
for (gene in tissuegenes){
  PlotBeforeAfter(array.exprs, array.adj, rna.seq.exprs)
}
dev.off()


