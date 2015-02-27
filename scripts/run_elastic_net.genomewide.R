# run_elastic_net.genomwide.R
# Jake Yeung
# February 26 2015

setwd("~/projects/tissue-specificity")

source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/RunElasticNet.R")

x <- read.table(file = "site_count_matrices//sitecounts_matrix_EPD", row.names = 1, header = TRUE)
conds <- c('Adr', 'Aorta', 'BFAT', 'Kidney', 'Liver', 'Lung', 'Mus')

for (cond in conds){
  y <- read.table(paste0("y_input_elastic_net/genome_wide//genome_wide.", cond, ".elasticinput"))
#   y <- read.table(paste0("y_input_elastic_net/7_conds_filtered26.", cond, ".elasticinput"))
#   x <- read.table("site_count_matrices/7_conds_filtered26.sitecounts")
  
  common.genes <- intersect(rownames(y), rownames(x))
  
  y <- as.matrix(y[common.genes, ])
  x <- as.matrix(x[common.genes, ])
  
  # Run elastic net ---------------------------------------------------------
  print(cond)
  start.time <- Sys.time()
  RunElasticNet(y, x, paste0("plots/elastic_net/genome_wide/7_conds_filtered26", cond), standardize = TRUE, alpha = 0.5)
  print(Sys.time() - start.time)
}
