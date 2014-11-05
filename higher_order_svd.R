# Jake Yeung
# November 4 2014
# higher_order_svd.R
#
# Explore data using higher order single value decomposition
# use rTensor package


# Library and custom functions --------------------------------------------

library(rTensor)

functions.dir <- 'scripts/functions'
source(file.path(functions.dir, 'hosvd.R'))  # modified from rTensor v 1.1 by James Li


# Define constants --------------------------------------------------------

N.TIMEPTS <- 24  # 24 time points, 2 hours per time point (over 48 hrs)
N.TISSUES <- 12  # 12 tissues

# Load data ---------------------------------------------------------------


# define dirs
data_dir <- "data"
fname <- "hogenesch_2014_rma.txt"    # data reprocessed by RMA package
# fname <- "hogenesch_2014_rma.ensemblnames.txt"

# load data
data_path <- file.path(data_dir, fname)
print(paste("Reading data from,", data_path, "May take a few a minutes."))
dat <- read.table(data_path)
print("Read data to memory.")

# Scale for PCA -----------------------------------------------------------

dat <- as.matrix(scale(dat))


# SVD ---------------------------------------------------------------------

# svd <- svd(dat)
# 
# plot(svd$v[, 3], type='o')

# Convert to Tensor -------------------------------------------------------

# 3-mode tensor.
# mode 1: genes (N ~ tens of thousands)
# mode 2: number of time points (N = 24)
# mode 3: number of tissues (N = 12)

tnsr <- fold(dat, rs=1, cs=c(2, 3), modes=c(nrow(dat), N.TIMEPTS, N.TISSUES))

# sanity check
# Below is matrix exprs of row 1 gene. Rows are time. Columns are tissues.
# tnsr[1, 1:N.TIMEPTS, 1:N.TISSUES]


# Higher Order SVD --------------------------------------------------------

# HOSVD, truncate after 100 because I get subscript out of bound errors?!

tnsr.hosvd <- hosvd(tnsr, rank=c(N.TIMEPTS * N.TISSUES, N.TIMEPTS, N.TISSUES))

# print dimensions of each orthogonal matrix, U
lapply(tnsr.hosvd$U, dim)


# Plot PCAs ---------------------------------------------------------------





