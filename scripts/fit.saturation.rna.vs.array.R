# fit.saturation.or.lm.R
# December 1 2014
# Jake Yeung

# Constants to change -----------------------------------------------------

pval <- 1.1  # for F-test
y.max <- 25  # for plotting
diagnostics.plot.out <- 'plots/diagnostics.saturation.fit.log2.probe.pdf'
before.after.plot.out <- 'plots/before.after.saturation.fit.log2.probe.pdf'
side.by.side.plot.out <- 'plots/rna.seq.saturation.fit.log2.probe.pdf'

mean.var.fit.outpath <- 'plots/noise.model.mean.var.pdf'

# Functions ---------------------------------------------------------------

# SATURATION FUNCTION

saturation <- function(params, r) params[2] + (params[1] * r) / (params[3] + r)
saturation.inv <- function(params, m) params[3] * (m - params[2]) / (params[2] + params[1] - m)

# LINEAR FUNCTIONS
linear <- function(params, x) params[1] + params[2] * x
linear.inv <- function(params, y) {
  if (is.na(params[2])){
    params[2] <- Inf
  }
  (y - params[1]) / params[2]
}

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names
source(file.path(functions.dir, 'RegressionFunctions.R'))  # optim function
source(file.path(functions.dir, 'GeneTissueCalculations.R'))
source(file.path(functions.dir, 'PlotFunctions.R'))

# Clock genes --------------------------------------------------------

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')
clockgenes <- c(clockgenes, "Sry")  # known to not be expressed via RNASeq


# Tissue genes ------------------------------------------------------------

tissuegenes <- c("Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
tissuegenes <- c(tissuegenes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
fname.rna.seq <- "rna_seq_deseq_counts_colnames_fixed.txt"
fname.array <- "array_exprs_colnames_fixed.txt"

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

# rna.seq.exprs <- log2(rna.seq.exprs + 1)

# log2 transform microarray -----------------------------------------------

array.exprs <- 2 ^ array.exprs

# Create subset and common genes so each microarray point has corresponding rnaseq
common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))
common.samples <- intersect(colnames(array.exprs), colnames(rna.seq.exprs))
array.exprs.subset.common.g <- array.exprs[common.genes, common.samples]  # log2 scale
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



# Fit saturation curve: 3 parameters --------------------------------------

# X-factor: how far from the min and max values do you allow?
# x <- 1: sensitive to noise
# large x: less accurate, but more conservative
x.factor <- 1

# init vals: linear model
slope0 <- 10
int0 <- 10

# init out list
fit.list <- vector(mode="list", length=length(common.genes))
names(fit.list) <- common.genes

for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  R.exprs <- unlist(rna.seq.exprs.common.g[gene, ])
  M.exprs <- unlist(array.exprs.subset.common.g[gene, ])
  M.full <- unlist(array.exprs[gene, ])
  
  # lower and upper bounds: saturation model
  bmin <- 0
  kmin <- 0
  amin <- max(M.full) * x.factor
  bmax <- min(M.full) / x.factor
  amax <- Inf
  kmax <- Inf
  # init vals: saturation model
  b0 <- 2^4
  k0 <- 2^7  # large to begin in linear regime
  a0 <- amin + 10  # expect it close to amin
  # lower and upper bounds: 
  slopemin <- 0
  intmin <- 0
  slopemax <- Inf
  intmax <- min(M.full) / 1.2
  # get variance as weights
  R.var <- predict(fit.noise, R.exprs)
  # adjust R.var so all values less than 0 take on smallest non-negative number
  R.var.min <- min(R.var[which(R.var > 0)])
  # if R.var.min is infinity, probably RNA-Seq has no expression, default weights to 1
  if (is.infinite(R.var.min)){
    warning(paste("Gene:", gene, "...adjusting infinites to 1 for variance."))
    R.var.min <- 1
  }
  R.var[which(R.var < 0)] <- R.var.min
  weights <- 1 / R.var
  
  fits <- tryCatch({
    
    fit.saturation <- nls(R.exprs ~ k * (M.exprs - b) / (b + a - M.exprs),
                         algorithm = "port",
                         start=list(a=a0, 
                                    b=b0,
                                    k=k0),
                         lower=list(a=amin,
                                    b=bmin,
                                    k=kmin),
                         upper=list(a=amax,
                                    b=bmax,
                                    k=kmax),
                         weights=weights)
    fits <- list(myfit=fit.saturation,
                 fit.used="saturation")
    
  }, error = function(e) {
    print(paste(gene, ', Error:', e))
    # fit.lm <- lm(M.exprs ~ R.exprs)
    fit.lm <- tryCatch({
      fit.lm <- nls(R.exprs ~ (M.exprs - int) / slope,
                    algorithm = "port",
                    start=list(int=int0,
                               slope=slope0),
                    lower=list(int=intmin,
                               slope=slopemin),
                    upper=list(int=intmax,
                               slope=slopemax),
                    weights=weights)
    }, error = function(e){
      warning(paste(gene, "Running normal lm. Need to convert."))
      # probably infinite slope, just fit normal lm
      fit.lm <- lm(R.exprs ~ M.exprs, 
                   weights=weights)
      return(fit.lm)
    })
    return(list(myfit=fit.lm,
                fit.used="lm"))
  })
  fit.list[gene] <- list(fits)
}

# Plot clock genes: diagnostics -------------------------------------------


pdf(diagnostics.plot.out)
par(mfrow = c(2,1))
for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  fit.used <- fit.list[[gene]][["fit.used"]]
  myfit <- fit.list[[gene]][["myfit"]]
  
  # Get vector of predicted values
  y <- GetFullR(gene, rna.seq.exprs, common.samples)
  x <- unlist(array.exprs[gene, ])
  x.predict <- seq(min(x)*0.9, max(x)*1.1, length.out=10*length(x))
  if (fit.used == "saturation"){
    y.hat <- saturation.inv(coef(myfit), x.predict)
  } else if (fit.used == "lm"){
    y.hat <- linear.inv(coef(myfit), x.predict)
  } else {
    warning("Neither saturation nor lm")
  }
  
  # plot 
  params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
  symbols <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=8, obs=1)
  sizes <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=0.25, obs=1)
  plot(x, y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       ylab="RNA-Seq DESeq-normalized counts",
       xlab="Microarray normal scale")
  lines(x.predict, y.hat)
  plot(log2(x + 1), log2(y + 1), main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       ylab="RNA-Seq DESeq-normalized counts (log2)",
       xlab="Microarray log2")
  lines(log2(x.predict + 1), log2(y.hat + 1))

}
dev.off()

# Adjust microarray -------------------------------------------------------

array.adj <- matrix(NA, nrow=nrow(array.exprs), ncol=ncol(array.exprs),
                    dimnames=list(rownames(array.exprs), 
                                  colnames(array.exprs)))

for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  fit.used <- fit.list[[gene]][["fit.used"]]  # either saturation or lm
  myfit <- fit.list[[gene]][["myfit"]]  
  
  array.exprs.gene <- as.matrix(array.exprs[gene, ])
  
  if (fit.used == "saturation"){
    array.adj.gene <- saturation.inv(coef(myfit), array.exprs.gene)
    
  } else if (fit.used == "lm"){
    array.adj.gene <- linear.inv(coef(myfit), array.exprs.gene)
  }
  array.adj[gene, ] <- array.adj.gene
}


# Plot before and after clockgenes ----------------------------------------

pdf(before.after.plot.out)
for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  tryCatch({
    PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, 
                    y.max=y.max, convert.log2=TRUE)
  }, warning=function(w){
    print(paste(gene, w))
    PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, 
                    y.max=y.max, convert.log2=TRUE)
  })
}
dev.off()


# Plot before and after against RNASeq ------------------------------------

pdf(side.by.side.plot.out)
for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  PlotAgainstRnaSeq(gene, log2(rna.seq.exprs.common.g + 1), log2(array.adj + 1), common.samples, y.max=y.max)
}
dev.off()
