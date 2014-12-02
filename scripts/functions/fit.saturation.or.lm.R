# fit.saturation.or.lm.R
# December 1 2014
# Jake Yeung

# Constants to change -----------------------------------------------------

pval <- 0.5  # for F-test
diagnostics.clock.plot.out <- 'plots/diagnostics.clock.lm.saturation.fit.log2.probe.pdf'
before.after.clock.plot.out <- 'plots/clockgenes.lm.saturation.fit.log2.probe.pdf'
diagnostics.tissue.plot.out <- 'plots/diagnostics.tissue.lm.saturation.fit.log2.probe.pdf'
before.after.tissue.plot.out <- 'plots/tissuegenes.lm.saturation.fit.log2.probe.pdf'

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



# Fit saturation curve: 3 parameters --------------------------------------

max.val <- Inf  # expected
# init vals
b0 <- 2
k0 <- 5  # large to begin in linear regime
a0 <- 15  # 13 is saturation i expect initially (a0 + b0)
# lower and upper bounds
bmin <- 0
kmin <- 0
amin <- 0
bmax <- 15
amax <- 1000
kmax <- 1000

# init out list
fit.list <- vector(mode="list", length=length(common.genes))
names(fit.list) <- common.genes

for (gene in c(clockgenes, tissuegenes)){
  R.exprs <- unlist(rna.seq.exprs.common.g[gene, ])
  M.exprs <- unlist(array.exprs.subset.common.g[gene, ])
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
                                    a=bmax,
                                    k=kmax))
    fit.lm <- lm(M.exprs ~ R.exprs)
    fits <- list(saturation=fit.saturation,
                 lm=fit.lm)
    
  }, error = function(e) {
    
    # print(paste('Error:', e))
    # error, try linear model
    fit.saturation <- NA
    fit.lm <- lm(M.exprs ~ R.exprs)
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
  x.predict <- seq(min(x)*0.8, max(x)*1.2, 0.1)
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
       xlab="RNA-Seq DESeq-normalized counts (log2)",
       ylab="Microarray log2")
  lines(x.predict, y.hat)
  plot(2^x - 1, 2^y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
       xlab="RNA-Seq DESeq-normalized counts",
       ylab="Microarray normal scale")
  lines(2^x.predict - 1, 2^y.hat)
}
dev.off()



