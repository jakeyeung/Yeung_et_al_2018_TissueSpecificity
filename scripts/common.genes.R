# common.genes.R

# define dirs
data.dir <- "data"
fname.array <- file.path(data.dir, "array_exprs_colnames_fixed.txt")
fname.array2 <- file.path(data.dir, "array_exprs_colnames_fixed2.txt")


# Load data ---------------------------------------------------------------

array1 <- read.table(fname.array)
array2 <- read.table(fname.array2)


# common.genes ------------------------------------------------------------

genes1 <- array1[, 1][2:nrow(array1)]
genes2 <- array2[, 1][2:nrow(array1)]

rm(array1)
rm(array2)


# Intersect ---------------------------------------------------------------

common.genes <- intersect(genes1, genes2)
