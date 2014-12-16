# Jake Yeung
# Dec 15 2014
# Project on time: but this time we do it INTEGRATED.

# Functions ---------------------------------------------------------------

scripts.dir <- "scripts"
funcs.dir <- file.path(scripts.dir, "functions")
source(file.path(funcs.dir, "DataHandlingFunctions.R"))
source(file.path(funcs.dir, "GetTissueTimes.R"))
source(file.path(funcs.dir, "RemoveProblemGenes.R"))
source(file.path(funcs.dir, "PlotFunctions.R"))  # for plotting complex functions
source(file.path(funcs.dir, "GetTopNValues.R"))

PlotVectorArrows <- function(complex.vector, exprs,
                             main="plot title", 
                             arrow.size=0.5,
                             dot.size=3){
  # Plot vector of complex numbers as additive arrows.
  # Please plot a blank plot first then run this.
  # e.g. plot(0, type="n", xlab="", ylab="", xlim=c(-max, max), ylim=c(-max, max), main=gene)
  # 
  # Args: 
  # complex.vector: vector of complex numbers whose sum means something...
  # exprs: complex number representing original expression
  #
  # Returns:
  # plot of sum of vector of complex numbers reaching its "destination" (i.e. sum of complex.vector)
  # plot "destination"
  points(exprs, cex=dot.size)
  abline(h=0)
  abline(v=0)
  
  # plot sum of arrows
  x0 <- 0  # start from origin
  y0 <- 0  # start from origin
  for (i in seq(1, length(complex.vector))){
    complex.val <- complex.vector[i]
    x1 <- x0 + Re(complex.val)
    y1 <- y0 + Im(complex.val)
    arrows(x0, y0, x1, y1, length=arrow.size)
    # redefine new starts as end of vector.
    x0 <- x1
    y0 <- y1
  }
}

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
array.normalized.fname <- "array.adj.0.07.txt"
array.normalized.path <- file.path(data.dir, array.normalized.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)


# Load file ---------------------------------------------------------------

array.normalized <- read.table(array.normalized.path)
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')



# Handle array ------------------------------------------------------------

array.normalized <- RemoveProblemGenes(array.normalized)

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns

# Transform to log2 scale -----------------------------------------------

array.normalized <- as.data.frame(log2(array.normalized + 1))
rna.seq.exprs <- as.data.frame(log2(rna.seq.exprs + 1))

# Define common genes -----------------------------------------------------

filtered.genes <- intersect(rownames(array.normalized), rownames(rna.seq.exprs))

array.normalized <- array.normalized[filtered.genes, ]
rna.seq.exprs <- rna.seq.exprs[filtered.genes, ]


# Merge array with RNASeq -------------------------------------------------

# merge to get 288 array samples, 96 rnaseq samples (in order)
merged.dat <- cbind(array.normalized, rna.seq.exprs)  

# Project frequency domain ------------------------------------------------

times.array <- as.integer(GetTimes(samp.names = colnames(array.normalized)))
times.rnaseq <- as.integer(GetTimes(samp.names = colnames(rna.seq.exprs)))

T <- 24
omega <- 2 * pi / T
n.tissues <- 12
n.timepts.array <- 24
interval.array <- 2
n.timepts.rnaseq <- 8
interval.rnaseq <- 6
tissues <- GetTissues(colnames(array.normalized))

# init output matrix
Y.gcs <- matrix(nrow = length(filtered.genes), ncol = length(tissues), dimnames = list(filtered.genes, tissues))

# Begin projection: iterate for each tissue.
for (tissue in tissues){
  # Get subset of merged.dat containing tissue
  merged.tissue <- as.matrix(merged.dat[, grepl(tissue, colnames(merged.dat))])
  # multiply by times.array and rna.seq array (unnormalizd fourier transform-ish)
  transformed <- merged.tissue %*% exp(-1i * omega * c(times.array, times.rnaseq))
  # normalize by number of data points
  transformed <- 1 / length(c(times.array, times.rnaseq)) * transformed
  Y.gcs[, tissue] <- transformed
}

Peek(Y.gcs)

# Let’s look at clock genes... --------------------------------------------

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')

for (gene in clockgenes){
  exprs <- Y.gcs[gene, ]
  max <- max(Mod(exprs))
  PlotComplex(exprs, 
              labels = tissues, 
              col = "HSV", 
              axis.min = -max, 
              axis.max = max, 
              main = gene, 
              rotate = pi / 2,
              add.text.plot = FALSE)
}


# Let’s SVD this ----------------------------------------------------------

s <- svd(Y.gcs)


# Rowname colname of SVD --------------------------------------------------

rownames(s$u) <- filtered.genes
rownames(s$v) <- tissues


# Add up outer products ---------------------------------------------------

# gene <- "Arntl"
# tissue <- "Liver"
n.components <- length(s$d)
n.components <- 3

outer.prods <- vector(length = n.components)
outer.prods <- list()

gene.list <- clockgenes
tissue.subset <- c("Liver", "Heart", "BS")

for (gene in gene.list){
  for (tissue in tissue.subset){
    outer.prod.vector <- vector(length=n.components)
    for (component in seq(1, n.components)){
      outer.prod <- s$d[component] * OuterComplex(s$u[gene, component], s$v[tissue, component]) 
      outer.prod.vector[component] <- outer.prod
    }
    outer.prods[[gene]][[tissue]] <- outer.prod.vector
  }  
}


# Plot as sum of vectors --------------------------------------------------

for (gene in gene.list){
  jmax <- 3 * sum(sapply(outer.prods[[gene]], function(x){
      max(Mod(x))
  }))
  plot(0, type="n", xlab="", ylab="", xlim=c(-jmax, jmax), ylim=c(-jmax, jmax), main=gene)
  for (tissue in tissue.subset){
    exprs <- Y.gcs[gene, tissue]
    PlotVectorArrows(outer.prods[[gene]][[tissue]], exprs, arrow.size = 0.1, dot.size=1.5)
  }
}


# Look at U and V matrices separately. ------------------------------------

components <- c(1, 2, 3, 4, 5)

# intersect.genes <- filtered.genes

for (component in components){
  top.genes <- GetTopNValues(Mod(s$u[, component]), N = 100)  # list of $vals $i
  # only take top.genes from U (makes it smaller, for visualization purposes)
  outer.prod.mat <- s$d[component] * OuterComplex(s$u[top.genes$i, component, drop = FALSE], t(s$v[, component]))
  # sort by phase angles
  # order.phase <- order(Arg(outer.prod.mat))
  # outer.prod.mat <- outer.prod.mat[order.phase, ]
  outer.prod.mat <- OrderPhaseMatrix(outer.prod.mat)
  PlotArgsMatrix(outer.prod.mat, main=paste("Component:", component))
  # intersect.genes <- intersect(intersect.genes, names(top.genes$vals))
  print(paste('-------------------Module', component))
  cat(paste0(names(top.genes$vals), collapse = '\n'))
  cat('\n')
}



