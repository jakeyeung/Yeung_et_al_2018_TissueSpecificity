# project.on.time.regulators.R
# Jake Yeung
# 18 Dec 2014
# SVD decomposition of matrices: on the level of regulators.

# Functions ---------------------------------------------------------------

scripts.dir <- "scripts"
funcs.dir <- file.path(scripts.dir, "functions")
source(file.path(funcs.dir, "DataHandlingFunctions.R"))
source(file.path(funcs.dir, "GetTissueTimes.R"))
source(file.path(funcs.dir, "RemoveProblemGenes.R"))
source(file.path(funcs.dir, 'FourierFunctions.R'))  # for Fourier stuff
source(file.path(funcs.dir, 'GetTopNValues.R'))
source(file.path(funcs.dir, "OuterComplex.R"))
source(file.path(funcs.dir, "OrderPhaseMatrix.R"))
source(file.path(funcs.dir, "PlotFunctions.R"))  # for plotting complex functions

get_TFs <- function(tf_vector){
  # Vector containing gene names (may be comma separated), get gene list
  # 
  # Args:
  # tf_vector: vector containing gene names (amy be comma separated)
  # with promoter (obtained from get_TFs_from_associations.py script)
  
  genes.list <- strsplit(tf_vector, split=",")
  
  # flatten list, return uniques
  genes.list <- unique(unlist(genes.list))
  
  return(genes.list)  
}

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
array.normalized.fname <- "array.adj.0.07.txt"
array.normalized.path <- file.path(data.dir, array.normalized.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)
tf.fname <- "motifs_and_TFs.list"
tf.path <- file.path(data.dir, tf.fname)


# Load file ---------------------------------------------------------------

array.normalized <- read.table(array.normalized.path)
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')
tf.mat <- read.table(tf.path, header=FALSE, row.names = 1, sep='\t')


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

# Define transcription factors --------------------------------------------

TFs <- get_TFs(as.vector(tf.mat[, 1]))
print(paste(length(TFs), "TFs found"))

tfs.filtered <- intersect(filtered.genes, TFs)
print(paste(length(tfs.filtered), "TFs found in RNA.Seq and array"))


# Filter merge data for TFs filtered --------------------------------------

merged.dat <- merged.dat[tfs.filtered, ]


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
tissues <- GetTissues(colnames(merged.dat))

# init output matrix
Y.gcs <- matrix(nrow = length(tfs.filtered), 
                ncol = length(tissues), 
                dimnames = list(tfs.filtered, tissues))

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

# Garbage -----------------------------------------------------------------

rm(rna.seq.exprs, array.normalized)


# SVD to see interesting components ---------------------------------------

s <- svd(Y.gcs)

# Rowname colname of SVD --------------------------------------------------

rownames(s$u) <- tfs.filtered
rownames(s$v) <- tissues

# Get heat map of components ----------------------------------------------


components <- c(1, 2, 3, 4, 5)

# intersect.genes <- filtered.genes

for (component in components){
  outfile <- paste0("results/module_", component, "_genes.txt")
  top.genes <- GetTopNValues(Mod(s$u[, component]), N = 100)  # list of $vals $i
  # only take top.genes from U (makes it smaller, for visualization purposes)
  outer.prod.mat <- s$d[component] * OuterComplex(s$u[top.genes$i, component, drop = FALSE], t(s$v[, component]))
  # sort by phase angles
  # order.phase <- order(Arg(outer.prod.mat))
  # outer.prod.mat <- outer.prod.mat[order.phase, ]
  outer.prod.mat <- OrderPhaseMatrix(outer.prod.mat)
  PlotArgsMatrix(outer.prod.mat, main=paste("Component:", component))
  # intersect.genes <- intersect(intersect.genes, names(top.genes$vals))
  #   print(paste('-------------------Module', component))
  #   cat(paste0(names(top.genes$vals), collapse = '\n'))
  #   cat('\n')
  # write all genes to file (ordered by values)
  top.genes.all <- GetTopNValues(Mod(s$u[, component]), N = length(s$u[, component]))
  # WriteGenesToFile(top.genes.all$vals, outfile)  # super SLOW damn R
}

# Look at U and V matrices separately -------------------------------------

pdf("plots/components.pdf")
for (component in components){
  print(paste('Plotting component', component))
  eigengene <- s$v[, component]
  max <- max(Mod(eigengene))
  PlotComplex(eigengene, 
              axis.min=-max, 
              axis.max=max,
              labels=tissues, 
              col = "HSV",
              main=paste("Tissue loadings for singular values:", component),
              rotate= pi / 2,
              add.text.plot = FALSE,
              jpch=20,
              threshold=0)
  
  eigensample <- s$u[, component]
  max <- max(Mod(eigensample))
  PlotComplex(eigensample, 
              axis.min=-max, 
              axis.max=max,
              labels=tfs.filtered, 
              col = "HSV",
              main=paste("Gene loadings for singular values:", component),
              rotate= pi / 2,
              add.text.plot = FALSE,
              jpch=1,
              threshold = 0.5 * max,
              verbose=TRUE)
}
dev.off()

gene.list <- c("Arntl", "Bhlhe40", "Clock", "Dbp", "Egr1", "Egr2", "Fosb", "Hlf", "Id1", "Junb", "Nfil3", "Npas2", "Srebf1", "Zbtb16")  # module 1
gene.list <- c("E2f2", "Egr1", "Ehf", "Ikzf1", "Myf6", "Myod1", "Myog", "Pax5", "Pax7", "Pitx2", "Snai3", "Spib", "Stat4", "Tcf7")  # module 2
gene.list <- c("Egr2", "Hnf4a", "Ikzf1", "Lef1", "Myf6", "Pax5", "Pax8", "Pitx2", "Pou2f2", "Spib", "Stat4")  # module 3
gene.list <- c("Atf5", "Egr1", "Fos", "Mafb", "Myc", "Npas2", "Rfx4", "Sox9", "Srebf1", "Zbtb16")  # module 5

for (gene in gene.list){
  gene.subset <- merged.dat[gene, 1:288]  # only array
  y.max <- max(gene.subset)
  plot(unlist(gene.subset), type='o', col=rep(1:length(tissues), each=24), main=gene, ylim=c(0, y.max))
  x.text <- 24 * seq(from = 0, length.out = length(tissues)) + 12
  y.text <- rep(0, length(tissues))
  text(x.text, y.text, labels = tissues)
  # plot also the complex value
  mod.max <- max(Mod(Y.gcs[gene, ]))
  PlotComplex(Y.gcs[gene, ], 
              labels = tissues, col = "HSV", 
              axis.min = -mod.max, 
              axis.max = mod.max, 
              main = gene, rotate = pi / 2, 
              add.text.plot = FALSE)
}



