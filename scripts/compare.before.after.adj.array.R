# variance_gain_from_adj.R
# Before and after adjusting microarray, do we get a gain in variance?
# Jake Yeung
# Dec 5 2014


# libraries ---------------------------------------------------------------

library(ggplot2)


# Constants ---------------------------------------------------------------

N.TISSUES <- 12
N.TIMEPTS <- 24
INTERVAL <- 2

# Functions ---------------------------------------------------------------

functions.dir <- file.path("scripts", "functions")
source(file.path(functions.dir, "RemoveProblemGenes.R"))
source(file.path(functions.dir, "GetTissueTimes.R"))
source(file.path(functions.dir, 'FourierFunctions.R'))  # for Fourier stuff


# Load files --------------------------------------------------------------

array.adj.path <- file.path(data.dir, "array.adj.0.07.txt")
data.path <- file.path(data.dir, "array_exprs_colnames_fixed.best.probe.selected.txt")
array.after <- read.table(array.adj.path)
array.before <- read.table(data.path, header=TRUE)  # needs handle rownames and colnames

# If needed handle duplicate row names ------------------------------------

rownames(array.before) <- make.names(array.before$gene, unique=TRUE)
# genes <- array.exprs$gene
# coords <- array.exprs$coordinates
# probeid <- array.exprs$ID_REF

drop.cols <- c("gene")
array.before <- array.before[, !(names(array.before) %in% drop.cols)]

# Remove problematic genes from array.after -------------------------------

array.after <- RemoveProblemGenes(array.after)

# Peek at dataframes ------------------------------------------------------

Peek(array.before)  # expect gene names as row names, tissues in columns
Peek(array.after)  # sample names are more readable...  


# log2 transform array.after ----------------------------------------------

array.after <- log2(array.after + 1)


# Get tissue and time -----------------------------------------------------

tissues <- GetTissues(samp.names = colnames(array.after))

# Project onto flat and rhythmic time components --------------------------

# RHYTHMIC
OMEGA <- 2 * pi / PERIOD
# OMEGA <- 0
array.before.time.proj <- ProjectToPeriodicTime(as.matrix(array.before), 
                                              N.TISSUES, 
                                              N.TIMEPTS, 
                                              INTERVAL, 
                                              OMEGA, 
                                              tissues)
array.after.time.proj <- ProjectToPeriodicTime(as.matrix(array.after),
                                               N.TISSUES,
                                               N.TIMEPTS,
                                               INTERVAL,
                                               OMEGA,
                                               tissues)


# Genes with highest expression in this projection? -----------------------

mean.exprs <- apply(Mod(array.before.time.proj), 1, mean)
var.exprs <- apply(Mod(array.before.time.proj), 1, var)
mean.exprs.after <- apply(Mod(array.after.time.proj), 1, mean)
var.exprs.after <- apply(Mod(array.after.time.proj), 1, var)


# List top  ---------------------------------------------------------------

(head(mean.exprs[order(mean.exprs, decreasing = TRUE)]))
(head(var.exprs[order(var.exprs, decreasing = TRUE)]))

(head(mean.exprs.after[order(mean.exprs.after, decreasing = TRUE)]))
(head(var.exprs.after[order(var.exprs.after, decreasing = TRUE)]))


# Do we get any variance gain? --------------------------------------------

(var.total.before <- apply(Mod(array.before.time.proj)^2, 2, sum))
(var.total.after <- apply(Mod(array.after.time.proj)^2, 2, sum))


# Plot data ---------------------------------------------------------------

var.df <- data.frame(tissue=c(names(var.total.before), names(var.total.after)), 
                     before.or.after=as.factor(rep(c("before", "after"), each = N.TISSUES)),
                     sum.Mod.Ygc1.sqr=c(unlist(var.total.before), unlist(var.total.after)))

ggplot(var.df, aes(x=tissue, y=sum.Mod.Ygc1.sqr, fill=before.or.after)) + geom_bar(stat="identity", position="dodge")
ggsave(filename = "plots/sum.Ygc.before.after.adj.pdf")