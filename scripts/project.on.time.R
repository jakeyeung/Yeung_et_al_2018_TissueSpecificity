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


# Transform to normal scale -----------------------------------------------

array.unadj <- as.data.frame((2^array.unadj))

# Handle array ------------------------------------------------------------

array.normalized <- RemoveProblemGenes(array.normalized)

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns


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

for (gene in top.genes.intersect){
  max <- max(Mod(array.unadj.time[gene, ]))
  PlotComplex(array.unadj.time, c(gene), 
              labels = tissues, 
              axis.min = -max, 
              axis.max = max, 
              main = gene, 
              add.text.plot = FALSE)
}

for (gene in top.genes.intersect){
  
  
  # lapply: dat is dat, name is experiment (for title)
  lapply(list(list(dat=array.unadj.time, name="array.unadj"),
              list(dat=array.normalized.time, name="array.normalized"),
              list(dat=rna.seq.exprs.time, name="rnaseq")),
         function(x, gene, labels){
           max <- max(Mod(x$dat[gene, ]))
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
         gene=c(gene),
         labels=tissues)
}

# array.unadj.time.mean <- sort(apply(array.unadj.time, 1, function(x){
#   return(median(Mod(x)))
# }), decreasing = TRUE)
# 
# (top.genes.unadj <- names(array.unadj.time.mean[1:20]))
# 
# for (gene in top.genes.unadj){
#   max <- max(Mod(array.unadj.time[gene, ]))
#   PlotComplex(array.unadj.time, c(gene), labels = tissues, 
#               axis.min = -max, axis.max = max, main = gene, add.text.plot = FALSE)
# }
# 
# array.adj.time.mean <- sort(apply(array.normalized.time, 1, function(x){
#   return(median(Mod(x)))
# }), decreasing = TRUE)
# 
# (top.genes.adj <- names(array.adj.time.mean[1:20]))
# 
# for (gene in top.genes.adj){
#   max <- max(Mod(array.unadj.time[gene, ]))
#   PlotComplex(array.unadj.time, c(gene), labels = tissues, 
#               axis.min = -max, axis.max = max, main = gene, add.text.plot = FALSE)
# }

