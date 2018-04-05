# run_elastic_net.R
# February 26 2015
# from Jingkui Wang

# install.packages("glmnet")
library(glmnet)
setwd("~/projects/tissue-specificity")
source("scripts/functions/RunElasticNet.R")

args <- commandArgs(trailingOnly = TRUE)
alpha = as.numeric(args[1])

# Load y and x ------------------------------------------------------------

conds <- c('Adr', 'Aorta', 'BFAT', 'Kidney', 'Liver', 'Lung', 'Mus')

for (cond in conds){
  y <- read.table(paste0("y_input_elastic_net/7_conds_filtered4.", cond, ".elasticinput"))
  x <- read.table("site_count_matrices/7_conds_filtered4.sitecounts")
  
  common.genes <- intersect(rownames(y), rownames(x))
  
  y <- as.matrix(y[common.genes, ])
  x <- as.matrix(x[common.genes, ])
  
  # Run elastic net ---------------------------------------------------------
  
  RunElasticNet(y, x, paste0("plots/elastic_net/7_conds_filtered4/7_conds_filtered4", cond), standardize = FALSE, alpha = alpha)
}
