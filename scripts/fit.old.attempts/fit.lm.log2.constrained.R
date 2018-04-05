# Jake Yeung
# Dec 3 2014
# Fit linear model: M on y axis, R on x-axis
# Add constrains

# Constants to change -----------------------------------------------------

y.max <- 25

# plot outputs
diagnostics.plot.out <- 'plots/diagnostics.lm.only.constrained.fit.log2.probe.pdf'
before.after.plot.out <- 'plots/before.after.lm.only.constrained.fit.log2.probe.pdf'
side.by.side.plot.out <- 'plots/side.by.side.lm.only.constrained.fit.log2.probe.against.rna.seq.pdf'

# table outputs
array.adj.output <- 'data/array.adjusted.lm.constrained.log2.txt'


# Functions ---------------------------------------------------------------

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
clockgenes <- c(clockgenes, "Sry")


# Tissue genes ------------------------------------------------------------

tissuegenes <- c("Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
tissuegenes <- c(tissuegenes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R


# Problematic genes -------------------------------------------------------

# has negative values when no constraints.
problematicgenes <- c("Tmem59l", "Saysd1", "H2.Q1", "Col25a1", "Clec18a", "X4930427A07Rik", 
                      "Ncan", "Crtam", "Fam43b", "Nphp4", "Nuf2", "Sox11", "Krt23", "Myo1h", 
                      "Mas1", "Cd207", "Tgif2", "Sdsl", "Gm8659", "Fsd1", "X2510049J12Rik", 
                      "X1600029I14Rik", "Syndig1l", "Cyp4x1", "E130309D14Rik", "Gm281", 
                      "Wdr27", "Daw1", "Tnfsf9", "Myo16")

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

# Remove genes that are not expressed across 96 samples -------------------

# min.exprs <- 5.5
# 
# expressed <- apply(rna.seq.exprs, 1, function(x){
#   if (max(x) < min.exprs){
#     return(FALSE)
#   } else {
#     return(TRUE)
#   }
# })
# rna.seq.exprs <- rna.seq.exprs[which(expressed), ]


# Create subsets ----------------------------------------------------------

# Create subset and common genes so each microarray point has corresponding rnaseq
common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))
common.samples <- intersect(colnames(array.exprs), colnames(rna.seq.exprs))
array.exprs.subset.common.g <- array.exprs[common.genes, common.samples]  # log2 scale
rna.seq.exprs.common.g <- rna.seq.exprs[common.genes, ]
array.common.g <- array.exprs[common.genes, ]
Peek(array.exprs.subset.common.g)
Peek(rna.seq.exprs.common.g)


# Fit lm model on log2 scale ----------------------------------------------

# init out list
fit.list <- vector(mode="list", length=length(common.genes))
names(fit.list) <- common.genes

# set init values
b0 <- 3  # background
a0 <- 0.5  # compression factor

for (gene in common.genes){
  R.exprs <- unlist(rna.seq.exprs.common.g[gene, ])  # 96 samples, matches RNASeq
  M.exprs <- unlist(array.exprs.subset.common.g[gene, ])  # 96 samples, matches RNASeq
  M.full <- unlist(array.exprs[gene, ])  # full 288 samples
  # set constraints
  bmin <- 0
  bmax <- min(M.full)
  amin <- 0
  amax <- Inf
  fits <- tryCatch({
    fit.lm <- nls(M.exprs ~ b + a * R.exprs,
                  algorithm = "port",
                  start=list(b=b0,
                             a=a0),
                  lower=list(b=bmin,
                             a=amin),
                  upper=list(b=bmax,
                             a=amax))
    fits <- list(lm=fit.lm)
    
  }, error = function(e) {
    warning(paste(gene, 'error:', e))
    # error, try linear model without constraints
    fit.lm <- lm(M.exprs ~ R.exprs)
    return(list(lm=fit.lm))
  })
  fit.list[gene] <- list(fits)
}

# Adjust microarray -------------------------------------------------------

array.adj <- matrix(NA, nrow=nrow(array.common.g), ncol=ncol(array.common.g),
                    dimnames=list(rownames(array.common.g), 
                                  colnames(array.common.g)))
for (gene in common.genes){
  myfit <- fit.list[[gene]][["lm"]]
  
  array.exprs.gene <- as.matrix(array.exprs[gene, ])
  array.adj.gene <- linear.inv(coef(myfit), array.exprs.gene)
  array.adj[gene, ] <- array.adj.gene
}

# Write microarray to output
write.table(array.adj, file=array.adj.output, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)

# Plot clock genes: diagnostics -------------------------------------------


pdf(diagnostics.plot.out)
par(mfrow = c(2,1))
for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  myfit <- fit.list[[gene]][["lm"]]
  
  # Get vector of predicted values
  x <- GetFullR(gene, rna.seq.exprs, common.samples)
  y <- unlist(array.exprs[gene, ])
  x.predict <- seq(min(x)*0.8, max(x)*1.2, 0.1)
  y.hat <- linear(coef(myfit), x.predict)
  
  # plot 
  params.str <- paste0(c('b=', ' a='), signif(as.vector(coef(myfit)), 2), collapse=",")
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


# Plot before and after clockgenes ----------------------------------------

pdf(before.after.plot.out)
for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  PlotBeforeAfter(gene, array.exprs, array.adj, rna.seq.exprs, y.max=y.max)
}
dev.off()


# Plot side by side -------------------------------------------------------

pdf(side.by.side.plot.out)
for (gene in c(clockgenes, tissuegenes, problematicgenes)){
  PlotAgainstRnaSeq(gene, rna.seq.exprs.common.g, array.adj, common.samples, y.max=y.max)
}
dev.off()

