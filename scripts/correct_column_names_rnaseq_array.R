# Correct column names between array and rna.seq

# define directories ------------------------------------------------------

# define dirs
data_dir <- "data"
fname.rna.seq <- "exprs_combined.genenames.txt"
fname.array <- "hogenesch_2014_rma.genenames.colnameordered.txt"


# Define outputs ----------------------------------------------------------

array.output <- file.path(data_dir, "array_exprs_colnames_fixed.txt")
rna.seq.output <- file.path(data_dir, "rna_seq_deseq_counts_colnames_fixed.txt")

# load data: RNASeq and microarray ----------------------------------------


# load data: rnaseq
rna.seq.path <- file.path(data_dir, fname.rna.seq)
print(paste("Reading data from,", rna.seq.path, "May take a few a minutes."))
rna.seq.exprs1 <- read.table(rna.seq.path, header=TRUE, sep='\t')
print("Read data to memory.")

# load data: microarray
array.path <- file.path(data_dir, fname.array)
print(paste("Reading data from,", array.path, "May take a few a minutes."))
array.exprs1 <- read.table(array.path, header=TRUE, sep='\t')  # has duplicate rownames
print("Read data to memory.")


# Get column names --------------------------------------------------------


array.cnames <- ShortenSampNames(colnames(array.exprs1)[2:ncol(array.exprs1)], show="tissue.time")

rna.cnames <- colnames(rna.seq.exprs1)[2:ncol(rna.seq.exprs1)]

t.names.array <- GetTissueNames(array.cnames, dat.type="array")

t.names.rnaseq <- GetTissueNames(rna.cnames, dat.type="rna.seq")


# Rename column names -----------------------------------------------------


# BEGIN: rename tissue names in rnaseq to match array
for (t.i in 1:length(t.names.array)){
  t.rnaseq <- t.names.rnaseq[t.i]
  t.array <- t.names.array[t.i]
  samp.names <- rna.cnames[grepl(t.rnaseq, rna.cnames)]
  # extract times for samp.names in tissue 
  samp.times <- GetTimes(samp.names, n.digits=2)
  # contactenate new tissue name to sample time
  samp.names.new <- paste0(t.array, samp.times)
  rna.cnames[grepl(t.rnaseq, rna.cnames)] <- samp.names.new
}
# END
print(rna.cnames)
print(array.cnames)


# Save to old column names ------------------------------------------------


colnames(array.exprs1)[2:ncol(array.exprs1)] <- array.cnames
colnames(rna.seq.exprs1)[2:ncol(rna.seq.exprs1)] <- rna.cnames

# match first element of array with rna
colnames(rna.seq.exprs1)[1] <- colnames(array.exprs1)[1]

# Write to file -----------------------------------------------------------


write.table(array.exprs1, file=array.output, quote=FALSE, sep="\t", row.names=FALSE)
write.table(rna.seq.exprs1, file=rna.seq.output, quote=FALSE, sep="\t", row.names=FALSE)
