# noise_model.R
# Jake Yeung
# November 26 2014
# Plot variance versus mean for each tissue and each gene.

# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names

# define directories ------------------------------------------------------

# define dirs
data_dir <- "data"
fname.rna.seq <- "rna_seq_deseq_counts_colnames_fixed.txt"


# Load file ---------------------------------------------------------------

# load data: rnaseq
rna.seq.path <- file.path(data_dir, fname.rna.seq)
print(paste("Reading data from,", rna.seq.path, "May take a few a minutes."))
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')
print("Read data to memory.")

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns


# Get column names --------------------------------------------------------

tissue.names <- GetTissueNames(colnames(rna.seq.exprs), dat.type="array")  # use array because these are "colnames_fixed" meaning they mtach array colnames

