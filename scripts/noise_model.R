# noise_model.R
# Jake Yeung
# November 26 2014
# Plot variance versus mean for each tissue and each gene.

library(ggplot2)

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

tissue.names <- GetTissueNames(colnames(rna.seq.exprs)) 


# Calculate mean and variance for each gene per tissue --------------------

N <- nrow(rna.seq.exprs) * length(tissue.names)  # one measurement for each gene for all tissues.
mean.var <- list(mean=matrix(NA, nrow=nrow(rna.seq.exprs), ncol=length(tissue.names),
                             dimnames = list(rownames(rna.seq.exprs), 
                                             tissue.names)), 
                 var=matrix(NA, nrow=nrow(rna.seq.exprs), ncol=length(tissue.names),
                            dimnames = list(rownames(rna.seq.exprs), tissue.names)))

for (j in 1:length(tissue.names)){
  tissue <- tissue.names[j]
  gene.tissue.exprs <- rna.seq.exprs[, grepl(tissue, colnames(rna.seq.exprs))]
  # calculate mean and var, by row
  exprs.mean <- apply(gene.tissue.exprs, 1, mean)
  exprs.var <- apply(gene.tissue.exprs, 1, var)
  # append to matrix 
  mean.var$mean[, j] <- exprs.mean
  mean.var$var[, j] <- exprs.var
}
str(mean.var)


# Make dataframe for ggplot2 ----------------------------------------------

mean.var.df <- data.frame(mean=as.vector(mean.var$mean), var=as.vector(mean.var$var))


# ggplot2 scatter ---------------------------------------------------------

ggplot(mean.var.df, aes(x=mean, y=var)) +
  geom_point(shape=46, alpha=0.5) + 
  scale_y_log10() +
  scale_x_log10()

# # Plot variance against mean ----------------------------------------------
# 
# set.seed(0)
# n.subset.rows <- ceiling(0.01 * nrow(rna.seq.exprs))
# subset.rows <- sample(rownames(rna.seq.exprs), n.subset.rows)
# x <- unlist(mean.var$mean[subset.rows, ])
# y <- unlist(mean.var$var[subset.rows, ])
# plot(x, y, log="xy", cex=0.1, pch='.')


# Fit loess ---------------------------------------------------------------

row.subset <- sample(nrow(mean.var.df), 0.01 * nrow(mean.var.df))
mean.var.df.subset <- mean.var.df[row.subset, ]

fit.loess <- loess(var ~ mean, mean.var.df.subset)
summary(fit.loess)
plot(mean.var.df.subset$mean, mean.var.df.subset$var, log="xy")
plot(fit.loess$x, predict(fit.loess), log="xy")
