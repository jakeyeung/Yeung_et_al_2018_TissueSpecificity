# run_elastic_net.R
# February 26 2015
# from Jingkui Wang

# install.packages("glmnet")
library(glmnet)
source("scripts/functions/RunElasticNet.R")

# Load y and x ------------------------------------------------------------

conds <- c('Adr', 'Aorta', 'BFAT', 'Kidney', 'Liver', 'Lung', 'Mus')

for (cond in conds){
  y <- read.table(paste0("y_input_elastic_net/7_conds_filtered26.", cond, ".elasticinput"))
  x <- read.table("site_count_matrices/7_conds_filtered26.sitecounts")
  
  common.genes <- intersect(rownames(y), rownames(x))
  
  y <- as.matrix(y[common.genes, ])
  x <- as.matrix(x[common.genes, ])
  
  # Run elastic net ---------------------------------------------------------
  
  RunElasticNet(y, x, paste0("plots/elastic_net/7_conds_filtered26_standardize", cond), standardize = FALSE, alpha = 0.5)
}
