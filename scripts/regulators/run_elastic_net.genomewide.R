# run_elastic_net.genomwide.R
# Jake Yeung
# February 26 2015

library(glmnet)

setwd("~/projects/tissue-specificity")

source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/RunElasticNet.R")

args <- commandArgs(trailingOnly = TRUE)
site_counts.fpath <- args[1]
output.dir <- args[2]
alpha <- as.numeric(args[3])

# site_counts.dir <- "site_count_matrices/GLM/SwissRegulon/N.MLE"
# output.dir <- "plots/elastic_net/genome_wide/GLM/SwissRegulon/N.MLE"
# alpha <- 0.05
# x <- read.table(file = file.path(site_counts.dir, "sitecount_matrix"), header = TRUE)
x <- read.table(file = site_counts.fpath, header = TRUE)
x.rownames <- make.names(x[, 1], unique = TRUE)
rownames(x) <- x.rownames
x$Gene.ID <- NULL

conds <- c('Liver')
for (cond in conds){
  y <- read.table(paste0("y_input_elastic_net/genome_wide//genome_wide.", cond, ".elasticinput"))
  head(rownames(x))
  head(rownames(y))
  common.genes <- intersect(rownames(y), rownames(x))
  print(paste0(length(common.genes), " genes in common."))
  y <- as.matrix(y[common.genes, ])
  x <- as.matrix(x[common.genes, ])
  
  # Run elastic net ---------------------------------------------------------
  print(cond)
  start.time <- Sys.time()
  RunElasticNet(y, x, file.path(output.dir, paste0("corrected_sitecounts_", cond)), standardize = TRUE, alpha = alpha)
  print(Sys.time() - start.time)
}

print(warnings())