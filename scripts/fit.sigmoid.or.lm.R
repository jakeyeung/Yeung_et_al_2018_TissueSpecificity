# Jake Yeung
# Dec 1 2014
# Fit linear model: M on y axis, R on x-axis

library(plyr)
library(ggplot2)


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


# Linear regression -------------------------------------------------------

coeff.mat <- matrix(NA, nrow=length(common.genes), ncol=2,
                    dimnames=list(common.genes, 
                                  c("slope", "intercept")))
for (gene in common.genes){
  R <- unlist(rna.seq.exprs.common.g[gene, ])
  M <- unlist(array.exprs.subset.common.g[gene, ])
  fit <- lm(M ~ R)  # R on y-axis, M on x-axis
  coeff.mat[gene, "intercept"] <- fit$coefficients[1]
  coeff.mat[gene, "slope"] <- fit$coefficients[2]
}

coeff.mat[clockgenes, ]

# Fit sigmoidal curve: general --------------------------------------------
# http://en.wikipedia.org/wiki/Generalised_logistic_function
# A: lower asymptote
# K: upper asymptote
# B: growth rate
# v > 0: affects near which asymptote max growth occurs
# Q: depends on the value Y(0)
# M: time to max growth if Q = v

# init vals
A0 <- 4
K0 <- 13
B0 <- 2
v0 <- 5
Q0 <- 10
M0 <- 7.5

# test
# gene <- "Elovl3"
# R <- unlist(rna.seq.exprs.common.g[gene, ])
# x <- unlist(array.exprs.subset.common.g[gene, ])

plot(R, sigmoid(A0, K0, B0, v0, Q0, M0, R))

gene.test <- c("Elovl3", "Npas2", "Arntl", "Per1")
for (gene in gene.test){
  R.exprs <- unlist(rna.seq.exprs.common.g[gene, ])
  M.exprs <- unlist(array.exprs.subset.common.g[gene, ])
  fit.sigmoid <- nls(M.exprs ~ A + ((K - A) / (1 + Q * exp(-B * (R.exprs - M))) ^ (1 / v)),
                     start=list(A=A0, 
                                K=K0,
                                B=B0,
                                v=v0,
                                Q=Q0,
                                M=M0))
  print(fit.sigmoid)
}

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
Bmax <- 6  # too high and the slope is super sensitive
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

for (gene in clockgenes){
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


# pdf('plots/diagnostics.lm.sigmoid.fit.log2.probe.pdf')
for (gene in clockgenes){
  # check if sigmoid fit was performed...
  fit.sigmoid <- fit.list[[gene]][["sigmoid"]]
  fit.lm <- fit.list[[gene]][["lm"]]
  
  # select either sigmoid or lm based on f.test or by default (i.e. sigmoid is NA)
  fit.select <- FTestSigmoidLinearModels(fit.lm, fit.sigmoid)
  fit.used <- fit.select$fit.used  # either sigmoid or lm
  myfit <- fit.select$myfit
  
  # Get vector of predicted values
  x <- unlist(rna.seq.exprs.common.g[gene, ])
  y <- unlist(array.exprs.subset.common.g[gene, ])
  x.predict <- seq(min(x), max(x), 0.1)
  if (fit.used == "sigmoid"){
    y.hat <- sigmoid(coef(myfit), x.predict)
  } else if (fit.used == "lm"){
    y.hat <- linear(coef(myfit), x.predict)
  } else {
    warning("Neither sigmoid nor lm")
  }
  
  # plot 
  params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
  plot(x, y, main=paste0("Gene=", gene, " Params=", params.str))
  lines(x.predict, y.hat)
  # plot(y, x, main=paste0("Gene=", gene, " Params=", params.str))
  # lines(y.hat, x.predict)
}
# dev.off()


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
  PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, y.max=16)
}
dev.off()


# Plot side by side: clock genes ------------------------------------------

pdf('plots/clock.genes.log2.lin.fit.by.probe.rna.seq.side.by.side.pdf')
for (gene in clockgenes){
  PlotAgainstRnaSeq(gene, rna.seq.exprs, array.adj, common.samples, y.max=16)  
}
dev.off()

# Plot tissue genes -------------------------------------------------------


pdf('plots/tissue.genes.log2.lin.fit.by.probe.pdf')
for (gene in tissuegenes){
  PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, y.max=16)
}
dev.off()


# Plot side by side: tissue genes -----------------------------------------

pdf('plots/tissue.genes.log2.lin.fit.by.probe.rna.seq.side.by.side.pdf')
for (gene in tissuegenes){
  PlotAgainstRnaSeq(gene, rna.seq.exprs, array.adj, common.samples, y.max=16)  
}
dev.off()

# Visualize tissue genes ---------------------------------------------------

# log(r) = a log(m) + b

pdf('plots/diagnostics.tissue.lm.fit.log2.probe.pdf')
for (gene in tissuegenes){
  slope <- coeff.mat[gene, "slope"]
  int <- coeff.mat[gene, "intercept"]
  PlotDiagnostics(gene, array.exprs, rna.seq.exprs, 
                  common.samples, slope, int)
}
dev.off()