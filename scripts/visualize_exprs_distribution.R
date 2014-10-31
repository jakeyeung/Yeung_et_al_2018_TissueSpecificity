# Jake Yeung
# October 31 2014
# visualize_exprs_distribution.R
# Plot distribution of gene expressions (may correlate with microarray?)


# define directories ------------------------------------------------------


# define dirs
data_dir <- "data"
fname.rna.seq <- "exprs_combined_rnaseq.txt"
fname.array <- "hogenesch_2014_rma.ensembl.noNA.genenames.txt"


# load data: RNASeq and microarray ----------------------------------------


# load data: rnaseq
rna.seq.path <- file.path(data_dir, fname.rna.seq)
print(paste("Reading data from,", rna.seq.path, "May take a few a minutes."))
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t', row.names=1)
print("Read data to memory.")

# load data: microarray
array.path <- file.path(data_dir, fname.array)
print(paste("Reading data from,", array.path, "May take a few a minutes."))
array.exprs <- read.table(array.path, header=TRUE, sep='\t')  # has duplicate rownames
print("Read data to memory.")


# Separate some columns in array data -------------------------------------

rownames(array.exprs) <- make.names(array.exprs$ensemblid, unique=TRUE)
# genes <- array.exprs$gene
# coords <- array.exprs$coordinates
# probeid <- array.exprs$ID_REF

drop.cols <- c("ensemblid", "gene", "coordinates", "ID_REF")
array.exprs <- array.exprs[, !(names(array.exprs) %in% drop.cols)]
str(array.exprs)

# Log2 Transform ----------------------------------------------------------

rna.seq.exprs <- log2(rna.seq.exprs)
array.exprs <- log2(array.exprs)

# Plot densities ----------------------------------------------------------

# plot density
plot(density(exprs.vector))

array.exprs.vector <- unlist(array.exprs)
plot(density(array.exprs.vector))
