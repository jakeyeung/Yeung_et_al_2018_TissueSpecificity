# fit.saturation.or.lm.R
# December 1 2014
# Jake Yeung

# Constants to change -----------------------------------------------------

pval <- 1e-5  # for F-test
diagnostics.clock.plot.out <- 'plots/diagnostics.clock.lm.constrained.saturation.fit.log2.probe.pdf'
before.after.clock.plot.out <- 'plots/clockgenes.lm.constrained.saturation.fit.log2.probe.pdf'
diagnostics.tissue.plot.out <- 'plots/diagnostics.tissue.lm.constrained.saturation.fit.log2.probe.pdf'
before.after.tissue.plot.out <- 'plots/tissuegenes.lm.constrained.saturation.fit.log2.probe.pdf'

# Functions ---------------------------------------------------------------

# SATURATION FUNCTION

saturation <- function(params, r) params[2] + (params[1] * r) / (params[3] + r)
saturation.inv <- function(params, m) params[3] * (m - params[2]) / (params[2] + params[1] - m)

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



# Fit saturation curve: 3 parameters --------------------------------------

# define init vals
max.val <- Inf  # expected
# init vals: saturation model
b0 <- 2^4
k0 <- 2^10  # large to begin in linear regime
a0 <- 2^13  # 13 is saturation i expect initially (a0 + b0)

# init vals: linear model
slope0 <- 10
int0 <- 10

# init out list
fit.list <- vector(mode="list", length=length(common.genes))
names(fit.list) <- common.genes

for (gene in c(clockgenes, tissuegenes)){
  R.exprs <- unlist(rna.seq.exprs.common.g[gene, ])
  M.exprs <- unlist(array.exprs.subset.common.g[gene, ])
  M.full <- unlist(array.exprs[gene, ])
  # lower and upper bounds: saturation model
  bmin <- 0
  kmin <- 0
  amin <- 0
  bmax <- min(M.full) - 1
  amax <- Inf
  kmax <- Inf
  # lower and upper bounds: 
  slopemin <- 0
  intmin <- 0
  slopemax <- Inf
  intmax <- min(M.full) - 1
  
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
                                    k=kmax))
    # fit.lm <- lm(M.exprs ~ R.exprs)
    fit.lm <- nls(M.exprs ~ int + slope * R.exprs,
                  algorithm = "port",
                  start=list(int=int0,
                             slope=slope0),
                  lower=list(int=intmin,
                             slope=slopemin),
                  upper=list(int=intmax,
                             slope=slopemax))
    
    fits <- list(saturation=fit.saturation,
                 lm=fit.lm)
    
  }, error = function(e) {
    
    # print(paste('Error:', e))
    # error, try linear model
    fit.saturation <- NA
    # fit.lm <- lm(M.exprs ~ R.exprs)
    fit.lm <- nls(M.exprs ~ int + slope * R.exprs,
                  algorithm = "port",
                  start=list(int=int0,
                             slope=slope0),
                  lower=list(int=intmin,
                             slope=slopemin),
                  upper=list(int=intmax,
                             slope=slopemax))
    return(list(saturation=fit.saturation, 
                lm=fit.lm))
    
  })
  fit.list[gene] <- list(fits)
}



# F-test on saturation and linear fit ----------------------------------------

fit.select.list <- vector(mode="list", length=length(common.genes))
names(fit.select.list) <- common.genes
for (gene in c(clockgenes, tissuegenes)){
  # check if saturation fit was performed...
  fit.saturation <- fit.list[[gene]][["saturation"]]
  fit.lm <- fit.list[[gene]][["lm"]]
  # select either saturation or lm based on f.test or by default (i.e. saturation is NA)
  fit.select <- FTestSigmoidLinearModels(fit.lm, fit.saturation, pval=pval, complex.model="saturation")
  fit.select.list[[gene]] <- fit.select
}

# Plot clock genes: diagnostics -------------------------------------------


pdf(diagnostics.clock.plot.out)
par(mfrow = c(2,1))
for (gene in clockgenes){
  fit.select <- fit.select.list[[gene]]
  fit.used <- fit.select$fit.used  # either saturation or lm
  myfit <- fit.select$myfit
  
  # Get vector of predicted values
  x <- GetFullR(rna.seq.exprs, common.samples)
  y <- unlist(array.exprs[gene, ])
  x.predict <- seq(min(x)*0.8, max(x)*1.2, 10)
  if (fit.used == "saturation"){
    y.hat <- saturation(coef(myfit), x.predict)
  } else if (fit.used == "lm"){
    y.hat <- linear(coef(myfit), x.predict)
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
  lines(x.predict, y.hat)
  plot(log2(x + 1), log2(y), main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts (log2)",
       ylab="Microarray log2")
  lines(log2(x.predict + 1), log2(y.hat))

}
dev.off()


# Plot tissue genes: diagnostics ------------------------------------------

pdf(diagnostics.tissue.plot.out)
par(mfrow = c(2,1))
for (gene in tissuegenes){
  fit.select <- fit.select.list[[gene]]
  fit.used <- fit.select$fit.used  # either saturation or lm
  myfit <- fit.select$myfit
  
  # Get vector of predicted values
  x <- GetFullR(rna.seq.exprs, common.samples)
  y <- unlist(array.exprs[gene, ])
  x.predict <- seq(min(x)*0.8, max(x)*1.2, 10)
  if (fit.used == "saturation"){
    y.hat <- saturation(coef(myfit), x.predict)
  } else if (fit.used == "lm"){
    y.hat <- linear(coef(myfit), x.predict)
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
  lines(x.predict, y.hat)
  plot(log2(x + 1), log2(y), main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts (log2)",
       ylab="Microarray log2")
  lines(log2(x.predict + 1), log2(y.hat))
  
}
dev.off()

# Adjust microarray -------------------------------------------------------

array.adj <- matrix(NA, nrow=nrow(array.exprs), ncol=ncol(array.exprs),
                    dimnames=list(rownames(array.exprs), 
                                  colnames(array.exprs)))

for (gene in c(clockgenes, tissuegenes)){
  fit.used <- fit.select.list[[gene]][["fit.used"]]
  myfit <- fit.select.list[[gene]][["myfit"]]  # either saturation or lm
  
  array.exprs.gene <- as.matrix(array.exprs[gene, ])
  
  if (fit.used == "saturation"){
    b <- coef(myfit)["b"]  # microarray cannot be less than b (bg)
    # adjust all microarray values less than A to A.
    array.exprs.gene[which(array.exprs.gene < b)] <- b
    array.adj.gene <- saturation.inv(coef(myfit), array.exprs.gene)
    
  } else if (fit.used == "lm"){
    int <- coef(myfit)[[1]]
    slope <- coef(myfit)[[2]]
    # adjust all microarray values less than int to int (bg)
    array.exprs.gene[which(array.exprs.gene < int)] <- int
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
  PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, 
                  y.max=y.max, convert.log2=TRUE)
}
dev.off()

# Plot before and after tissuegenes ----------------------------------------

pdf(before.after.tissue.plot.out)
for (gene in tissuegenes){
  PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, 
                  y.max=y.max, convert.log2=TRUE)
}
dev.off()

