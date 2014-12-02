# Jake Yeung
# Dec 1 2014
# Fit linear model: M on y axis, R on x-axis

library(plyr)
library(ggplot2)


# Constants to change -----------------------------------------------------

pval <- 1e-100  # for F-test
y.max <- 20
diagnostics.clock.plot.out <- 'plots/diagnostics.clock.lm.only.fit.log2.probe.pdf'
before.after.clock.plot.out <- 'plots/clockgenes.lm.only.fit.log2.probe.pdf'
side.by.side.clock.plot.out <- 'plots/clockgenes.lm.only.fit.log2.probe.against.rna.seq.pdf'
diagnostics.tissue.plot.out <- 'plots/diagnostics.tissue.lm.only.fit.log2.probe.pdf'
before.after.tissue.plot.out <- 'plots/tissuegenes.lm.only.fit.log2.probe.pdf'
side.by.side.tissue.plot.out <- 'plots/tissuegenes.lm.only.fit.log2.probe.against.rna.seq.pdf'


# Functions ---------------------------------------------------------------

# SIGMOID FUNCTIONS 
# http://psg.hitachi-solutions.com/masterplex/blog/the-4-parameter-logistic-4pl-nonlinear-regression-model
sigmoid <- function(params, x) ((params[1] - params[4]) / (1 + ((x / params[3]) ^ params[2]))) + params[4]
sigmoid.inv <- function(params, y) params[3] * ((params[1] - y) / (y - params[4])) ^ (1 / params[2])
#  A is the MFI (Mean Fluorescent Intensity)/RLU (Relative Light Unit) value for the minimum asymptote
#  B is the Hill slope
#  C is the concentration at the inflection point
#  D is the MFI/RLU value for the maximum asymptote
#  E is the asymmetry factor

# LINEAR FUNCTIONS
linear <- function(params, x) params[1] + params[2] * x
linear.inv <- function(params, y) (y - params[1]) / params[2]

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
clockgenes <- c(clockgenes, "Sry")


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

rna.seq.exprs <- log2(rna.seq.exprs + 1)
# Create subset and common genes so each microarray point has corresponding rnaseq
common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))
common.samples <- intersect(colnames(array.exprs), colnames(rna.seq.exprs))
array.exprs.subset.common.g <- array.exprs[common.genes, common.samples]  # log2 scale
rna.seq.exprs.common.g <- rna.seq.exprs[common.genes, ]
Peek(array.exprs.subset.common.g)
Peek(rna.seq.exprs.common.g)


# Fit sigmoidal curve: 4 params --------------------------------------------

# init vals
A0 <- 4
B0 <- 1
C0 <- 8
D0 <- 13
Amin <- 2
Bmin <- 0
Cmin <- 3
Dmin <- 13
Amax <- Dmin
Bmax <- 5  # too high and the slope is super sensitive
Cmax <- 100
Dmax <- Cmax

# test
# gene <- "Per1"
# R <- unlist(rna.seq.exprs.common.g[gene, ])
# x <- unlist(array.exprs.subset.common.g[gene, ])
# 
# plot(R, sigmoid.4(A0, B0, C0, D0, R), main='4 param.')
# plot(R, sigmoid.5(A0, B0, C0, D0, R), main='5 param: E=1')
# plot(R, x)

# init out list
fit.list <- vector(mode="list", length=length(common.genes))
names(fit.list) <- common.genes

for (gene in c(clockgenes, tissuegenes)){
  R.exprs <- unlist(rna.seq.exprs.common.g[gene, ])
  M.exprs <- unlist(array.exprs.subset.common.g[gene, ])
  fits <- tryCatch({
    
    fit.sigmoid <- nls(M.exprs ~ ((A - D) / (1 + ((R.exprs / C) ^ B))) + D,
                       algorithm = "port",
                       start=list(A=A0, 
                                  B=B0,
                                  C=C0,
                                  D=D0),
                       lower=list(A=Amin,
                                  B=Bmin,
                                  C=Cmin,
                                  D=Dmin),
                       upper=list(A=Amax,
                                  B=Bmax,
                                  C=Cmax,
                                  D=Dmax))
    fit.lm <- lm(M.exprs ~ R.exprs)
    fits <- list(sigmoid=fit.sigmoid,
                 lm=fit.lm)
    
  }, error = function(e) {
    
    # print(paste('Error:', e))
    # error, try linear model
    fit.sigmoid <- NA
    fit.lm <- lm(M.exprs ~ R.exprs)
    return(list(sigmoid=fit.sigmoid, 
                lm=fit.lm))
    
  })
  fit.list[gene] <- list(fits)
}


# F-test on sigmoid and linear fit ----------------------------------------

fit.select.list <- vector(mode="list", length=length(common.genes))
names(fit.select.list) <- common.genes
for (gene in c(clockgenes, tissuegenes)){
  # check if sigmoid fit was performed...
  fit.sigmoid <- fit.list[[gene]][["sigmoid"]]
  fit.lm <- fit.list[[gene]][["lm"]]
  # select either sigmoid or lm based on f.test or by default (i.e. sigmoid is NA)
  fit.select <- FTestSigmoidLinearModels(fit.lm, fit.sigmoid, pval=pval)
  fit.select.list[[gene]] <- fit.select
}


# Plot clock genes: diagnostics -------------------------------------------


pdf(diagnostics.clock.plot.out)
par(mfrow = c(2,1))
for (gene in clockgenes){
  fit.select <- fit.select.list[[gene]]
  fit.used <- fit.select$fit.used  # either sigmoid or lm
  myfit <- fit.select$myfit
  
  # Get vector of predicted values
  x <- GetFullR(rna.seq.exprs, common.samples)
  y <- unlist(array.exprs[gene, ])
  x.predict <- seq(min(x)*0.8, max(x)*1.2, 0.1)
  if (fit.used == "sigmoid"){
    y.hat <- sigmoid(coef(myfit), x.predict)
  } else if (fit.used == "lm"){
    y.hat <- linear(coef(myfit), x.predict)
  } else {
    warning("Neither sigmoid nor lm")
  }
  
  # plot 
  params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
  symbols <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=8, obs=1)
  sizes <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=0.25, obs=1)
  plot(x, y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts (log2)",
       ylab="Microarray log2")
  lines(x.predict, y.hat)
  plot(2^x - 1, 2^y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts",
       ylab="Microarray normal scale")
  lines(2^x.predict - 1, 2^y.hat)
}
dev.off()


# Adjust microarray -------------------------------------------------------

array.adj <- matrix(NA, nrow=nrow(array.exprs), ncol=ncol(array.exprs),
                    dimnames=list(rownames(array.exprs), 
                                  colnames(array.exprs)))
for (gene in c(clockgenes, tissuegenes)){
  fit.used <- fit.select.list[[gene]][["fit.used"]]
  myfit <- fit.select.list[[gene]][["myfit"]]  # either sigmoid or lm
  
  array.exprs.gene <- as.matrix(array.exprs[gene, ])
  
  if (fit.used == "sigmoid"){
    A <- coef(myfit)["A"]  # microarray cannot be less than A (bg)
    # adjust all microarray values less than A to A.
    array.exprs.gene[which(array.exprs.gene < A)] <- A
    array.adj.gene <- sigmoid.inv(coef(myfit), array.exprs.gene)

  } else if (fit.used == "lm"){
    int <- coef(myfit)[[1]]
    slope <- coef(myfit)[[2]]
    if (is.na(slope)){
      # if NA, assume it is infinity
      slope <- Inf
    }
    array.adj.gene <- (array.exprs.gene - int) / slope
  }
  array.adj[gene, ] <- array.adj.gene
}


# Plot before and after clockgenes ----------------------------------------

pdf(before.after.clock.plot.out)
for (gene in clockgenes){
  PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, y.max=y.max)
}
dev.off()


# Plot side by side -------------------------------------------------------

pdf(side.by.side.clock.plot.out)
for (gene in clockgenes){
  PlotAgainstRnaSeq(gene, rna.seq.exprs.common.g, array.adj, common.samples, y.max=y.max)
}
dev.off()

# Plot tissue genes: diagnostics -------------------------------------------


pdf(diagnostics.tissue.plot.out)
par(mfrow = c(2,1))
for (gene in tissuegenes){
  fit.select <- fit.select.list[[gene]]
  fit.used <- fit.select$fit.used  # either sigmoid or lm
  myfit <- fit.select$myfit
  
  # Get vector of predicted values
  x <- GetFullR(rna.seq.exprs, common.samples)
  y <- unlist(array.exprs[gene, ])
  x.predict <- seq(min(x)*0.8, max(x)*1.2, 0.1)
  if (fit.used == "sigmoid"){
    y.hat <- sigmoid(coef(myfit), x.predict)
  } else if (fit.used == "lm"){
    y.hat <- linear(coef(myfit), x.predict)
  } else {
    warning("Neither sigmoid nor lm")
  }
  
  # plot 
  params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
  symbols <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=8, obs=1)
  sizes <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=0.25, obs=1)
  plot(x, y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts (log2)",
       ylab="Microarray log2")
  lines(x.predict, y.hat)
  plot(2^x - 1, 2^y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts",
       ylab="Microarray normal scale")
  lines(2^x.predict - 1, 2^y.hat)
}
dev.off()

# Plot before and after tissuegenes ----------------------------------------

pdf(before.after.tissue.plot.out)
for (gene in tissuegenes){
  PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, y.max=y.max)
}
dev.off()


# Plot side by side tissue genes ------------------------------------------

pdf(side.by.side.tissue.plot.out)
for (gene in tissuegenes){
  PlotAgainstRnaSeq(gene, rna.seq.exprs.common.g, array.adj, common.samples, y.max=y.max)
}
dev.off()
