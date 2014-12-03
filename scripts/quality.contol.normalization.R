# Jake Yeung
# quality.control.normalization.R
# Dec 2 2014

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
normalized.array.fname <- "array.adj.saturationfit.conservative.txt"
normalized.array.path <- file.path(data.dir, normalized.array.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)

# Load file ---------------------------------------------------------------

normalized.array <- read.table(normalized.array.path)
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')


# Handle rnaseq -----------------------------------------------------------

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns

# Remove rows with NA -----------------------------------------------------

normalized.array <- na.omit(normalized.array)


# How many have negative values? ------------------------------------------

negs <- apply(normalized.array, 1, function(x){
  if (min(x) < 0){
    return(1)
  } else {
    return(0)
  }
})

problem.genes <- names(negs[which(negs == 1)])


# Remove problem genes ----------------------------------------------------

print(paste0("Removing " length(problem.genes), " problem genes."))

normalized.array <- normalized.array[which(!rownames(normalized.array) %in% problem.genes), ]

# Plot density in log2 ----------------------------------------------------

par(mfrow = c(2, 1))
plot(density(unlist(log2(normalized.array + 1))))
plot(density(unlist(log2(rna.seq.exprs[common.genes, ] + 1))))
par(mfrow = c(1, 1))