# Jake Yeung
# November 25 2014
# fit curve



# Functions ---------------------------------------------------------------


functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'DataHandlingFunctions.R'))  # for peeking at Data


# define directories ------------------------------------------------------

# define dirs
data_dir <- "data"

fname.array <- file.path(data_dir, "array_exprs_colnames_fixed.txt")
fname.rna.seq <- file.path(data_dir, "rna_seq_deseq_counts_colnames_fixed.txt")


# Load data ---------------------------------------------------------------


# load data: rnaseq
rna.seq.path <- file.path(fname.rna.seq)
print(paste("Reading data from,", rna.seq.path, "May take a few a minutes."))
rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')
print("Read data to memory.")

# load data: microarray
array.path <- file.path(fname.array)
print(paste("Reading data from,", array.path, "May take a few a minutes."))
array.exprs <- read.table(array.path, header=TRUE, sep='\t')  # has duplicate rownames
print("Read data to memory.")


# Handle duplicate rownames: RNASEQ ------------------------------------

rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]

Peek(rna.seq.exprs)  # expect gene names as row names, tissues in columns

# Handle duplicate rownames: MICROARRAY -------------------------------------

# first column contained gene names
rownames(array.exprs) <- make.names(array.exprs$gene, unique=TRUE)

drop.cols <- c("gene")
array.exprs <- array.exprs[, !(names(array.exprs) %in% drop.cols)]
Peek(array.exprs)


# log2 transform RNASEQ ---------------------------------------------------

rna.seq.exprs <- log2(rna.seq.exprs + 1)


# Get intersection of genes -----------------------------------------------

common.genes <- intersect(rownames(array.subset), rownames(rna.seq.exprs))
array.common.g <- as.matrix(array.exprs[common.genes, ])
rna.seq.common.g <- as.matrix(rna.seq.exprs[common.genes, ])

# Subset exprs ------------------------------------------------------------

array.subset <- array.common.g[, colnames(rna.seq.exprs)]


# Create function R in function of M --------------------------------------

a.init <- 2^0
b.init <- 2^4

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12', 
                'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')

for (gene in clockgenes){
  R <- rna.seq.common.g[gene, ]
  M <- array.subset[gene, ]
  # x[2] is b
  # x[1] is a
  S <- function(x) sum((R - log(((2 ^ M - x[2]) / x[1]) + 1)) ^ 2)
  min.exprs <- min((2 ^ M))
  print(gene)
  print(S(c(a.init, b.init)))
  a.b <- optim(c(a.init, b.init), S, method="L-BFGS-B",
               lower=c(2^-2, 0),
               upper=c(2^2, min.exprs))
  print(a.b)
  
  
  # Plot results ------------------------------------------------------------
  
  f <- function(x) log(((2 ^ x - a.b$par[2]) / a.b$par[1]) + 1)
  x <- seq(from=0, to=12, length.out=100)
  y <- f(x)
  plot(M, R, main=paste(gene, 'a=', signif(a.b$par[1], 2), 'b=', signif(a.b$par[2], 2), 'converge=', a.b$convergence), 
       xlab="log2 microarray (M)", ylab="log2 rnaseq (R)")
  lines(x, y)
}

M.Cry1 <- array.subset['Cry1', ]
R.Cry1 <- rna.seq.common.g['Cry1', ]
x.space <- seq(from=0, to=12, length.out=100)
f <- function(x) log(((2 ^ x.space - x[2]) / x[1]) + 1)
y.Cry1 <- f(c(0.1, 0))
plot(M.Cry1, R.Cry1)
lines(x.space, y.Cry1)

