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
  max <- max(Mod(exprs))
  PlotComplex(Y.gcs, gene.list = c(gene), 
              labels = tissues, col = "HSV", 
              axis.min = -max, 
              axis.max = max, 
              main = gene, 
              add.text.plot = FALSE)
}
