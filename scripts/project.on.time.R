# project.on.time.R
# do Fourier on different datasets
# Jake Yeung
# Dec 11 2014

# Functions ---------------------------------------------------------------

scripts.dir <- "scripts"
funcs.dir <- file.path(scripts.dir, "functions")
source(file.path(funcs.dir, "DataHandlingFunctions.R"))
source(file.path(funcs.dir, "GetTissueTimes.R"))
source(file.path(funcs.dir, "RemoveProblemGenes.R"))
source(file.path(funcs.dir, "FourierFunctions.R"))
source(file.path(funcs.dir, "GetTissueTimes.R"))
source(file.path(funcs.dir, "PlotFunctions.R"))  # for plotting complex functions

OuterComplex <- function(u, v){
  # Given two complex vectors, get the outer product.
  # 
  # given u = x + i*y and v = x2 + i*y2, we want to get
  # u * Conj(v) = (x + i*y) * (x2 - i*y2) = x*x2 - y*y2 - 1i (x*y2 - x2*y)
  # 
  # Note in R: using outer(u, v) returns (x + i*y) * (x2 + i*y2) which
  # is not equivalent.
  # 
  # Args:
  # u -> complex number/vector
  # v -> complex number/vector
  # 
  # Returns
  # o -> outer product of u and v (u * Conj(v))
  # 
  o <- u %*% Conj(v)
  return(o)
}

GetMaxAcrossSets <- function(...){
  datasets <- list(...)
  max.across.sets <- max(sapply(datasets, function(x){
    max(Mod(x))
  }))
}

GetTopGenes <- function(dat, N="all"){
  # Args:
  # 
  # dat: complex matrix (projected to time)
  # N: number of genes to return. If N="all", returns all
  
  dat.mean <- sort(apply(dat, 1, function(x){
    return(median(Mod(x)))
  }), decreasing = TRUE)
  
  if (N=="all"){
    return(top.genes <- names(dat.mean))
  } else
    return(top.genes <- names(dat.mean[1:N]) )
}

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
array.unadj.fname <- "array_exprs_colnames_fixed.best.probe.selected.txt"
array.unadj.path <- file.path(data.dir, array.unadj.fname)
array.normalized.fname <- "array.adj.0.07.txt"
array.normalized.path <- file.path(data.dir, array.normalized.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)


# Load file ---------------------------------------------------------------

array.unadj <- read.table(array.unadj.path, header=TRUE, sep='\t')
array.normalized <- read.table(array.normalized.path)
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')



# Handle array unadjusted -------------------------------------------------

# first column contained gene names
rownames(array.unadj) <- make.names(array.unadj$gene, unique=TRUE)

drop.cols <- c("gene")
array.unadj <- array.unadj[, !(names(array.unadj) %in% drop.cols)]
Peek(array.unadj)


# Handle array ------------------------------------------------------------

array.normalized <- RemoveProblemGenes(array.normalized)

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns

# Transform to log2 scale -----------------------------------------------

# array.unadj <- as.data.frame((2^array.unadj))
array.normalized <- as.data.frame(log2(array.normalized + 1))
rna.seq.exprs <- as.data.frame(log2(rna.seq.exprs + 1))

# Define common genes -----------------------------------------------------

filtered.genes <- intersect(rownames(array.normalized), rownames(rna.seq.exprs))

array.normalized <- array.normalized[filtered.genes, ]
rna.seq.exprs <- rna.seq.exprs[filtered.genes, ]
array.unadj <- array.unadj[filtered.genes, ]


# Do SVD ------------------------------------------------------------------

T <- 24
omega <- 2 * pi / T
n.tissues <- 12
n.timepts <- 24
interval <- 2
n.timepts.rnaseq <- 8
interval.rnaseq <- 6
tissues <- GetTissues(colnames(array.normalized))

array.normalized.time <- ProjectToPeriodicTime(as.matrix(array.normalized), 
                                               n.tissues, 
                                               n.timepts, 
                                               interval, 
                                               omega, 
                                               tissues)

array.unadj.time <- ProjectToPeriodicTime(as.matrix(array.unadj), 
                                          n.tissues, 
                                          n.timepts, 
                                          interval, 
                                          omega, 
                                          tissues)

rna.seq.exprs.time <- ProjectToPeriodicTime(as.matrix(rna.seq.exprs),
                                            n.tissues, 
                                            n.timepts.rnaseq, 
                                            interval.rnaseq, 
                                            omega, 
                                            tissues)


# What does this look like? -----------------------------------------------

# Plot top genes by average.
top.genes.unadj <- GetTopGenes(array.unadj.time, N=100)
top.genes.normalized <- GetTopGenes(array.normalized.time, N=100)
top.genes.rnaseq <- GetTopGenes(rna.seq.exprs.time, N=100)

top.genes.intersect <- intersect(top.genes.unadj, top.genes.normalized)
top.genes.intersect <- intersect(top.genes.intersect, top.genes.rnaseq)


# Plot genes: separate PDF for each dataset -------------------------------


# Plot top genes: separate PDF for each data set.
lapply(list(list(dat=array.unadj.time, name="array.unadj"),
            list(dat=array.normalized.time, name="array.normalized"),
            list(dat=rna.seq.exprs.time, name="rnaseq")),
       function(x, gene.list, labels, 
                array.unadj.time, array.normalized.time, rna.seq.exprs.time){
         pdf(file.path("plots", paste0(x$name, ".pdf")))
         for (gene in gene.list){
           array.unadj.exprs <- array.unadj.time[gene, ]
           array.normalized.exprs <- array.normalized.time[gene, ]
           rna.seq.exprs <- rna.seq.exprs.time[gene, ]
           max <- GetMaxAcrossSets(array.unadj.exprs, 
                                   array.normalized.exprs, 
                                   rna.seq.exprs)
           min <- -max
           PlotComplex(x$dat,
                       gene,
                       labels,
                       col="HSV",
                       min,
                       max,
                       main=paste(gene, x$name),
                       add.text.plot=FALSE) 
         }
         dev.off()
       }, 
       gene.list=top.genes.intersect,
       labels=tissues,
       array.unadj.time=array.unadj.time,
       array.normalized.time=array.normalized.time,
       rna.seq.exprs.time=rna.seq.exprs.time)


# Plot genes: ONE pdf for all 3 datasets ----------------------------------


# Plot top genes: ONE pdf (subplots) for all 3 datasets.
pdf(file.path("plots", "top.genes.complex.pdf"))
par(mfrow = c(2,2))
for (gene in top.genes.intersect){
  lapply(list(list(dat=array.unadj.time, name="array.unadj"),
              list(dat=array.normalized.time, name="array.normalized"),
              list(dat=rna.seq.exprs.time, name="rnaseq")),
         function(x, gene, labels, 
                  array.unadj.time, array.normalized.time, rna.seq.exprs.time){
           array.unadj.exprs <- array.unadj.time[gene, ]
           array.normalized.exprs <- array.normalized.time[gene, ]
           rna.seq.exprs <- rna.seq.exprs.time[gene, ]
           max <- GetMaxAcrossSets(array.unadj.exprs, 
                                   array.normalized.exprs, 
                                   rna.seq.exprs)
           # add a factor to max so it's not on the "edge"
           max <- 1.2 * max
           min <- -max
           
           PlotComplex(x$dat,
                       gene,
                       labels,
                       col="HSV",
                       min,
                       max,
                       main=paste(gene, x$name),
                       add.text.plot=FALSE) 
         },
         gene = gene,
         labels = tissues,
         array.unadj.time = array.unadj.time,
         array.normalized.time = array.normalized.time,
         rna.seq.exprs.time = rna.seq.exprs.time)
  frame()  # plot blank to fill in the 2x2 matrix of subplots.
}
dev.off()



# Do SVD to see how the vector breaks down... -----------------------------

s.unadj <- svd(array.unadj.time)
s.norm <- svd(array.normalized.time)
s.rnaseq <- svd(rna.seq.exprs.time)
s.norm.orig <- svd(as.matrix(array.normalized))

# Plot my screeplots
lapply(list(list(dat=s.unadj, name="array.unadj"),
            list(dat=s.norm, name="array.normalized"),
            list(dat=s.rnaseq, name="rnaseq")), 
       function(lst){
         x <- lst$dat
         main <- paste("Screeplot", lst$name)
         plot(x$d / sum(x$d), type='o', main=main, 
              xlab="Singular values",
              ylab=expression(D[ii] / sum(D_ii)))
})


# Add rownames and colnames to my SVDs ------------------------------------

rownames(s.unadj$u)  <- rownames(array.unadj.time)
rownames(s.norm$u) <- rownames(array.normalized.time)
rownames(s.rnaseq$u) <- rownames(rna.seq.exprs.time)
rownames(s.norm.orig$u) <- rownames(array.normalized)

rownames(s.unadj$v) <- colnames(array.unadj.time)
rownames(s.norm$v) <- colnames(array.normalized.time)
rownames(s.rnaseq$v) <- colnames(rna.seq.exprs.time)
rownames(s.norm.orig$v) <- colnames(array.normalized)



# Let’s look at clock genes... --------------------------------------------

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')

gene <- c("Nr1d1")
# gene <- rownames(array.normalized.time)
tissue <- "Liver"
component <- 1
s <- s.norm
orig <- array.normalized.time
n.components <- ncol(orig)

(exprs.complex <- orig[gene, tissue])

exprs.uv.vector <- vector(length = n.components)

for (component in seq(1, n.components)){
  # exprs.uv <- s$d[component] * outer(s$u[gene, component], s$v[tissue, component])
  # exprs.uv <- s$u[gene, component] * s$d[component] * Conj(s$v[tissue, component])
  exprs.uv <- s$d[component] * OuterComplex(s$u[gene, component], s$v[tissue, component])
  exprs.uv.vector[component] <- exprs.uv
}
print(sum(exprs.uv.vector))






