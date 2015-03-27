# fastICA_rhythmic_genes.R
# Jake Yeung
# Date: 2015-03-20


# Functions ---------------------------------------------------------------

library(fastICA)
library(JADE)
library(doMC)
source('scripts/functions/LoadArrayRnaSeq.R')
source('scripts/functions/DataHandlingFunctions.R')
source('scripts/functions/ProjectToFrequency.R')
source('scripts/functions/SvdFunctions.R')
source('scripts/functions/ConvertLongToWide.R')
# Main --------------------------------------------------------------------

dat <- LoadArrayRnaSeq()

dat.split <- split(dat, dat$tissue)


# Remove WFAT from analysis -----------------------------------------------

# WFAT seems to be show a lot of epidydmal-specific genes, let's remove 

# it from analysis to prevent strangely oscillating genes showing up
dat.split$WFAT <- NULL


# Project my data ---------------------------------------------------------

omega <- 2 * pi / 24
omegas <- GetOmegas()

start.time <- Sys.time()
if (getDoParWorkers() == 1){
  registerDoMC(40)
}
print(paste("Parallel processing on:", getDoParWorkers(), "cores"))
dat.split.proj <- lapply(dat.split, function(x){
  ddply(x, .(gene), ProjectToFrequency, 
        my.omega = omega, 
        normalize = FALSE, 
        rhythmic.only = FALSE, 
        pval.cutoff = 1, 
        .parallel = TRUE)
})
print(Sys.time() - start.time)


# Add tissue information into each list -----------------------------------

for (tissue in names(dat.split.proj)){
  dat.split.proj[[tissue]]$tissue <- tissue
}

# Combine data ------------------------------------------------------------

dat.proj <- do.call(rbind, dat.split.proj)

# long to wide conversion
dat.wide <- ConvertLongToWide(dat.proj, measurement.var = "exprs.transformed")

# Complete cases only. This removes NaN rows.
dat.wide <- dat.wide[complete.cases(dat.wide), ]

# Run fastICA -------------------------------------------------------------

n.comp <- 20
ic <- fastICA(dat.wide, n.comp, alg.typ = "parallel",
              fun = c("logcosh","exp"), alpha = 1.0, method = "C",
              row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE,
              w.init = NULL)
ic.cjd <- cjd(as.matrix(dat.wide), eps = 1e-06, maxiter = 100)

ic.jade <- JADE(as.numeric(dat.wide), n.comp = 11, eps = 1e-06, maxiter = 100, na.action = na.fail)
