# fit.saturation.or.lm.R
# December 1 2014
# Jake Yeung

# Constants to change -----------------------------------------------------

pval <- 1.1  # for F-test between saturation model and linear model
y.max <- 25  # for plotting
diagnostics.plot.out <- 'plots/diagnostics.saturation.slope.adj.0.7.RERUN.pdf'
before.after.plot.out <- 'plots/before.after.saturation.slope.adj.0.7.RERUN.pdf'
side.by.side.plot.out <- 'plots/overlap.with.rnaseq.saturation.slope.adj.0.7.RERUN.pdf'

mean.var.fit.outpath <- 'plots/mean.var.probes.pdf'

array.adj.out <- 'data/array.adj.saturationfit.conservative.slope.adj.0.7.RERUN.txt'
fit.results.out <- 'data/saturation.fit.results.slope.adj.0.7.RERUN.RData'
fit.select.results.out <- 'data/saturation.fit.select.results.slope.adj.0.7.RERUN.RData'

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
source(file.path(functions.dir, 'FitSaturationCurveFindProblemGenes.R'))

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


# Problem genes -----------------------------------------------------------

# interesting cases
problematicgenes <- c("Tmem59l", "Saysd1", "H2.Q1", "Col25a1", "Clec18a", "X4930427A07Rik", 
                      "Ncan", "Crtam", "Fam43b", "Nphp4", "Nuf2", "Sox11", "Krt23", "Myo1h", 
                      "Mas1", "Cd207", "Tgif2", "Sdsl", "Gm8659", "Fsd1", "X2510049J12Rik", 
                      "X1600029I14Rik", "Syndig1l", "Cyp4x1", "E130309D14Rik", "Gm281", 
                      "Wdr27", "Daw1", "Tnfsf9", "Myo16")

# problem genes from array data without probe selection
# problematicgenes <- c(problematicgenes, 'Gm4793', 'Sumo1', 'Prima1', 'Cinp', 'Rnf180', 'Gja3', 'Tnf', 
#                   'Plin4', 'Rpl9', 'Dapl1', 'Tsc22d2', 'Clspn', 'Gnat3', 'Nupl2', 
#                   'Tmem213', 'Ssbp1', 'Zfp746', 'Abhd13', 'Gtpbp3', 'Adam5', 
#                   'Lonrf1', 'Gm16380', 'Kri1', 'Csnk2a1.ps', 'Dusp9')

# problem genes from array data with probe selection
problematicgenes <- c("Abhd13", "Adam5", "Cinp", "Clspn", "Csnk2a1.ps", 
                      "Dapl1", "Dusp9", "Gja3", "Gm16380", "Gm4793", "Gnat3", "Gtpbp3", 
                      "Kri1", "Lonrf1", "Nupl2", "Plin4", "Prima1", "Rnf180", "Rpl9", 
                      "Ssbp1", "Tmem213", "Tnf", "Tsc22d2", "Zfp746")

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
fname.rna.seq <- "rna_seq_deseq_counts_colnames_fixed.txt"
fname.array <- "array_exprs_colnames_fixed.best.probe.selected.txt"

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

mean.var.df <- GetMeanVarByTissues(exprs=array.exprs, tissue.names)

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
x <- seq(min(array.exprs), max(array.exprs), 100)  # plot full range
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


# Crazy loop to find best starting slope. ---------------------------------

slopes <- seq(0.06, by=0.001, length.out=50)
slope.n.neg.genes <- data.frame(n.neg.genes=rep(NA, length(slopes)), row.names=slopes)
for (slope0 in slopes){
  fit.results.out.slope <- paste0('data/sat.fit.', slope0, '.RData')
  fit.select.results.out.slope <- paste0('data/sat.fit.select.', slope0, '.RData')
  array.adj.out.slope <- paste0('data/array.adj.', slope0, '.txt')
  diagnostics.plot.out.slope <- paste0('plots/diagnostics.plot.', slope0, '.pdf')
  n.neg.genes <- FitSaturationCurveFindProblemGenes(slope0, 
                                                    fit.results.out.slope,
                                                    fit.select.results.out.slope,
                                                    array.adj.out.slope,
                                                    diagnostics.plot.out.slope)
  print(paste(n.neg.genes, 'negative genes.'))
  slope.n.neg.genes[toString(slope0), 1] <- n.neg.genes
}



# Fit saturation curve: 3 parameters --------------------------------------

# x.factor: make fit not so close to saturation points
# x.factor = 1: best fit
# x.factor > 1: force saturation points farhter from data points
x.factor <- 1.2
# bg.factor: make fit not so close to background
bg.factor <- 0.8

# init vals: linear model
slope0 <- 0.07
int0 <- 10

# init out list
fit.list <- vector(mode="list", length=length(common.genes))
names(fit.list) <- common.genes

for (gene in common.genes){
  R.exprs <- unlist(rna.seq.exprs.common.g[gene, ])
  M.exprs <- unlist(array.exprs.subset.common.g[gene, ])
  M.full <- unlist(array.exprs[gene, ])
  
  # lower and upper bounds: saturation model
  bmin <- 0
  kmin <- 0
  amin <- max(M.full) * x.factor
  bmax <- min(M.full) * bg.factor
  amax <- Inf
  kmax <- Inf
  # init vals: saturation model
  b0 <- 2^4
  k0 <- 2^6  # large to begin in linear regime
  a0 <- amin + 10  # expect it close to amin
  # lower and upper bounds: 
  slopemin <- 0
  intmin <- 0
  slopemax <- Inf
  intmax <- min(M.full) * bg.factor
  # get variance as weights
  M.var <- predict(fit.noise, M.exprs)
  # adjust M.var so all values less than 0 take on smallest non-negative number
  M.var.min <- min(M.var[which(M.var > 0)])
  # if M.var.min is infinity, probably RNA-Seq has no expression, default weights to 1
  if (is.infinite(M.var.min)){
    warning(paste("Gene:", gene, "...adjusting infinites to 1 for variance."))
    M.var.min <- 1
  }
  M.var[which(M.var < 0)] <- M.var.min
  weights <- 1 / M.var
  
  fits <- tryCatch({
    
    fit.saturation <- nls(M.exprs ~ b + (a * R.exprs) / (k + R.exprs),
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
    fit.lm <- FitLmConstraint(M.exprs, R.exprs, weights, 
                              int0, slope0, 
                              intmin, slopemin, 
                              intmax, slopemax)
    
    fits <- list(saturation=fit.saturation,
                 lm=fit.lm)
    
  }, error = function(e) {
    fit.saturation <- NA
    # fit.lm <- lm(M.exprs ~ R.exprs)
    fit.lm <- FitLmConstraint(M.exprs, R.exprs, weights, 
                              int0, slope0, 
                              intmin, slopemin, 
                              intmax, slopemax)
    
    return(list(saturation=fit.saturation, 
                lm=fit.lm))
    
  })
  fit.list[gene] <- list(fits)
}

# Save fit results to .RData ----------------------------------------------

save(fit.list, file = fit.results.out)

# F-test on saturation and linear fit ----------------------------------------

fit.select.list <- vector(mode="list", length=length(common.genes))
names(fit.select.list) <- common.genes
for (gene in common.genes){
  # check if saturation fit was performed...
  fit.saturation <- fit.list[[gene]][["saturation"]]
  fit.lm <- fit.list[[gene]][["lm"]]
  # select either saturation or lm based on f.test or by default (i.e. saturation is NA)
  fit.select <- FTestSigmoidLinearModels(fit.lm, fit.saturation, pval=pval, complex.model="saturation")
  fit.select.list[[gene]] <- fit.select
}

save(fit.select.list, file = fit.select.results.out)

# Plot clock genes: diagnostics -------------------------------------------

pdf(diagnostics.plot.out)
par(mfrow = c(2,1))
for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  fit.select <- fit.select.list[[gene]]
  fit.used <- fit.select$fit.used  # either saturation or lm
  myfit <- fit.select$myfit
  
  # Get vector of predicted values
  x <- GetFullR(gene, rna.seq.exprs, common.samples)
  y <- unlist(array.exprs[gene, ])
  # x.predict <- seq(min(x)*0.8, max(x)*1.2, length.out=10*length(x))
  y.predict <- seq(min(y)*0.8, max(y), length.out=10*length(y))
  if (fit.used == "saturation"){
    x.hat <- saturation.inv(coef(myfit), y.predict)
  } else if (fit.used == "lm"){
    x.hat <- linear.inv(coef(myfit), y.predict)
  } else {
    warning("Neither saturation nor lm")
  }
  
  # plot 
  params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
  symbols <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=8, obs=1)
  sizes <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=0.25, obs=1)
  plot(x, y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts",
       ylab="Microarray normal scale")
  lines(x.hat, y.predict)
  plot(log2(x + 1), log2(y), main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts (log2)",
       ylab="Microarray log2")
  lines(log2(x.hat + 1), log2(y.predict))
  
}
dev.off()


# Adjust microarray -------------------------------------------------------

array.adj <- matrix(NA, nrow=nrow(array.exprs), ncol=ncol(array.exprs),
                    dimnames=list(rownames(array.exprs), 
                                  colnames(array.exprs)))

for (gene in common.genes){
  fit.used <- fit.select.list[[gene]][["fit.used"]]
  myfit <- fit.select.list[[gene]][["myfit"]]  # either saturation or lm
  
  array.exprs.gene <- as.matrix(array.exprs[gene, ])
  
  if (fit.used == "saturation"){
    array.adj.gene <- saturation.inv(coef(myfit), array.exprs.gene)
    
  } else if (fit.used == "lm"){
    array.adj.gene <- linear.inv(coef(myfit), array.exprs.gene)
  }
  array.adj[gene, ] <- array.adj.gene
}


# Write to file -----------------------------------------------------------

write.table(array.adj, file = array.adj.out, 
            quote=FALSE, sep='\t',
            row.names=TRUE, col.names=NA)


# Check for problematic genes ---------------------------------------------

# How many have negative values? ------------------------------------------

negs <- apply(array.adj, 1, function(x){
  if (min(x) < 0){
    return(1)
  } else {
    return(0)
  }
})

problem.genes <- names(negs[which(negs == 1)])

print(paste('slope:', slope0, 'has', length(problem.genes), 'problem genes.'))

# Plot diagnostics for problem genes --------------------------------------

pdf(diagnostics.plot.out)
par(mfrow = c(2,1))
for (gene in c(problem.genes)){
  fit.select <- fit.select.list[[gene]]
  fit.used <- fit.select$fit.used  # either saturation or lm
  myfit <- fit.select$myfit
  
  # Get vector of predicted values
  x <- GetFullR(gene, rna.seq.exprs, common.samples)
  y <- unlist(array.exprs[gene, ])
  # x.predict <- seq(min(x)*0.8, max(x)*1.2, length.out=10*length(x))
  y.predict <- seq(min(y)*0.8, max(y), length.out=10*length(y))
  if (fit.used == "saturation"){
    x.hat <- saturation.inv(coef(myfit), y.predict)
  } else if (fit.used == "lm"){
    x.hat <- linear.inv(coef(myfit), y.predict)
  } else {
    warning("Neither saturation nor lm")
  }
  
  # plot 
  params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
  symbols <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=8, obs=1)
  sizes <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=0.25, obs=1)
  plot(x, y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts",
       ylab="Microarray normal scale")
  lines(x.hat, y.predict)
  plot(log2(x + 1), log2(y), main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts (log2)",
       ylab="Microarray log2")
  lines(log2(x.hat + 1), log2(y.predict))
  
}
dev.off()


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

# Garbage collect ---------------------------------------------------------

rm(fit.select.list)

# Garbage collection ------------------------------------------------------

rm(fit.list)
