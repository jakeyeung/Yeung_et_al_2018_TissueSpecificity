# Jake Yeung
# October 31 2014
# visualize_exprs_distribution.R
# Plot distribution of gene expressions (may correlate with microarray?)

library(ggplot2)
library(reshape2)

# Constants ---------------------------------------------------------------

epsilon <- 1  # RNA-Seq log2 transform to prevent infinities
N.RNASEQ.TIMEPTS <- 8
N.RNASEQ.INTERVAL <- 6
N.ARRAY.TIMEPTS <- 24
N.ARRAY.INTERVAL <- 2

# Define output directory (USER CHANGEABLE) -------------------------------

plot_dir <- "plots"
pdf_file <- "batch.bgcorrected.normalized.original.pdf"
pdf_path <- file.path(plot_dir, pdf_file)


# Functions ---------------------------------------------------------------

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data
source(file.path(functions.dir, 'GetTissueNames.R'))
source(file.path(functions.dir, 'SampleNameHandler.R'))  # make sample names

# define directories ------------------------------------------------------

# define dirs
data_dir <- "data"
fname.rna.seq <- "exprs_combined.genenames.txt"
fname.array <- "hogenesch_2014_rma.genenames.txt"
# fname.array <- "all.hogenesch.seq.background.normalized.genenames.cut.filter.txt"
# fname.array <- "all.hogenesch.background.nonormalize.genenames.txt"  # normalize MANUALLY by setting NORMALIZE to TRUE
NORMALIZE <- FALSE


# load data: RNASeq and microarray ----------------------------------------


# load data: rnaseq
rna.seq.path <- file.path(data_dir, fname.rna.seq)
print(paste("Reading data from,", rna.seq.path, "May take a few a minutes."))
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')
print("Read data to memory.")

# load data: microarray
array.path <- file.path(data_dir, fname.array)
print(paste("Reading data from,", array.path, "May take a few a minutes."))
array.exprs <- read.table(array.path, header=TRUE, sep='\t')  # has duplicate rownames
print("Read data to memory.")
Peek(array.exprs)

# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$Gene.ID, unique=TRUE)

drop.cols <- c("Gene.ID")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns

# Log2 Transform ----------------------------------------------------------

rna.seq.exprs <- log2(rna.seq.exprs + epsilon)
# array exprs already in log2

# Handle duplicate rownames: MICROARRAY -------------------------------------

# first column contained gene names
rownames(array.exprs) <- make.names(array.exprs$gene, unique=TRUE)
# genes <- array.exprs$gene
# coords <- array.exprs$coordinates
# probeid <- array.exprs$ID_REF

drop.cols <- c("gene")
array.exprs <- array.exprs[, !(names(array.exprs) %in% drop.cols)]
Peek(array.exprs)

# Normalize microarray ----------------------------------------------------

if (NORMALIZE == TRUE){
  library(preprocessCore)
  array.exprs <- normalize.quantiles(as.matrix(array.exprs), copy=FALSE)
}


# Handle array exprs tissue names -----------------------------------------

colnames(array.exprs) <- ShortenSampNames(colnames(array.exprs), show="tissue.time")
Peek(array.exprs)

pdf(file=pdf_path)

# Plot densities ----------------------------------------------------------
par(mfrow=c(1, 2))  # turn on side by side plot to visualize array and seq
# plot density
plot(density(unlist(rna.seq.exprs)))
plot(density(unlist(array.exprs)))
# turn off side by side pot
par(mfrow=c(1, 1))


# Scatterplot -------------------------------------------------------------

# subset rows so they have common rows
common.genes <- intersect(rownames(array.exprs), rownames(rna.seq.exprs))

# colnames to grep. Grep only times from RNA-Seq (last two digits of colnames)
# Times in RNA-Seq are subset of times in array.exprs
rna.seq.times <- sapply(colnames(rna.seq.exprs[1 : N.RNASEQ.TIMEPTS]), 
                        function(x){
  x <- substr(x, nchar(x) - 1, nchar(x))
  return(x)
})  # 8 time pts
# add asterisks in front to as our grep pattern, grabs everything that ends with RNA-seq times
rna.seq.times <- paste0('*', rna.seq.times)

array.exprs.subset <- array.exprs[, grepl(paste0(rna.seq.times, collapse='|'), colnames(array.exprs))]

# plot scatter
set.seed(0)
FRACTION.SAMPLES <- 0.05
random.sample <- sample(common.genes, FRACTION.SAMPLES * length(common.genes))
# plot(unlist(rna.seq.exprs[random.sample, ]), unlist(array.exprs.subset[random.sample, ]))

# use ggplot2, so we melt first then plot. We care about common genes only.
rna.melt <- melt(as.matrix(rna.seq.exprs[random.sample, ]))
array.melt <- melt(as.matrix(array.exprs.subset[random.sample, ]))
# check gene names are same
# all(rna.melt$Var1 == array.melt$Var1)  # TRUE

# combine rna and array into one df
rna.array.long <- data.frame("gene" = rna.melt$Var1, 
                             "sample" = rna.melt$Var2, 
                             "rna.exprs" = rna.melt$value, 
                             "array.exprs" = array.melt$value)
# make DF longer, melt rna and array exprs into a single column.
rna.array.long <- melt(rna.array.long)
colnames(rna.array.long) <- c("gene", "sample", "rnaseq.or.array", "exprs")
head(rna.array.long)


m <- ggplot(rna.array.long, aes(y=exprs[which(rnaseq.or.array == "rna.exprs")],
                            x=exprs[which(rnaseq.or.array == "array.exprs")]))
m <- m + geom_point(shape=1, alpha=0.25)    # Use hollow circles
print(m)  # to file


# Histogram and density ---------------------------------------------------


# Add histogram of samples: RNA-Seq and Array
ggplot(rna.array.long, aes(x=exprs, fill=rnaseq.or.array)) + 
  geom_histogram(binwidth=0.5) +
  facet_wrap(~ sample)

# Density
ggplot(rna.array.long, aes(x=exprs, fill=rnaseq.or.array)) + 
  geom_density(alpha=0.5) +
  facet_wrap(~ sample)


# Scatterplot “replicates" ------------------------------------------------


# Scatterplot: by "24 hour cycles"
# create list of "replicate pairs", i.e. samples separated by 24 hrs
period.pairs <- list()
for (i in 1 : (length(rna.seq.times) / 2 - 1)) {
  period.pairs[[i]] <- list(rep1=rna.seq.times[i], rep2=rna.seq.times[i + 4]) 
}

# plot rep1 with rep2: both RNA.Seq
for (i in 1 : (length(rna.seq.times) / 2 - 1)) {
  m <- ggplot(rna.array.long, aes(x=exprs[which(rnaseq.or.array == "rna.exprs" & grepl(period.pairs[[i]]$rep1, sample))],
                             y=exprs[which(rnaseq.or.array == "rna.exprs" & grepl(period.pairs[[i]]$rep2, sample))]))
  m <-  m + geom_point(shape=1, alpha=0.05)    # Use hollow circles
  m <- m + xlab(paste("Time:", period.pairs[[i]]$rep1, "hr"))
  m <- m + ylab(paste("Time:", period.pairs[[i]]$rep2, "hr"))
  m <- m + labs(title="RNA-Seq samples between 24 hours")
  print(m)
}

# plot rep1 with rep2: both microarray
for (i in 1 : (length(rna.seq.times) / 2 - 1)) {
  m <- ggplot(rna.array.long, aes(x=exprs[which(rnaseq.or.array == "array.exprs" & grepl(period.pairs[[i]]$rep1, sample))],
                             y=exprs[which(rnaseq.or.array == "array.exprs" & grepl(period.pairs[[i]]$rep2, sample))]))
  m <- m + geom_point(shape=1, alpha=0.05)    # Use hollow circles
  m <- m + xlab(paste("Time:", period.pairs[[i]]$rep1, "hr"))
  m <- m + ylab(paste("Time:", period.pairs[[i]]$rep2, "hr"))
  m <- m + labs(title="Array between 24 hours")
  print(m)
}

# plot rep1 with rep1: RNA Seq vs Microarray
for (t in rna.seq.times) {
  m <- ggplot(rna.array.long, aes(x=exprs[which(rnaseq.or.array == "rna.exprs" & grepl(t, sample))],
                                  y=exprs[which(rnaseq.or.array == "array.exprs" & grepl(t, sample))]))
  m <-  m + geom_point(shape=1, alpha=0.05)    # Use hollow circles
  m <- m + xlab(paste("RNA Exprs. Time:", t, "hr"))
  m <- m + ylab(paste("Array Exprs. Time:", t, "hr"))
  m <- m + labs(title=paste("RNA Exprs and Array exprs at", t, "hr"))
  print(m)
}


# Boxplots ----------------------------------------------------------------

rna.array.long.filtered <- subset(rna.array.long)
m <- ggplot(rna.array.long.filtered, aes(x=sample, y=exprs))
m <- m + geom_boxplot() + facet_wrap(~rnaseq.or.array)
print(m)

rna.array.long.filtered <- subset(rna.array.long, rnaseq.or.array="array.exprs")
m <- ggplot(rna.array.long.filtered, aes(x=sample, y=exprs))
m <- m + geom_boxplot()
print(m)

dev.off()

rm(rna.melt)
rm(array.melt)

print(paste("Plots saved in:", pdf_path))

