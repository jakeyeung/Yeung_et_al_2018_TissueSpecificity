# 2015-10-11
# Jake Yeung
# Loop through gene.Robjs containing
# topN models after running nconds2.
# rbind everything to create a large matrix. Save to file.

library(data.table)

start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")
source("scripts/functions/AppendListFunctions.R")


# Load --------------------------------------------------------------------

gene.dir <- "/home/yeung/projects/nconds_results/cat_fits"
outpath <- "/home/yeung/projects/nconds_results/fits_long_11_tiss_max_3.bug_fixed.rbindlist.Robj"

fits.list <- expandingList()
for (gene.robj in list.files(gene.dir)){
  load(file.path(gene.dir, gene.robj))  # fits.long.gene
  fits.list$add(fits.long.gene)
}

fits.long.lst <- fits.list$as.list()

fits.long <- rbindlist(l = fits.long.lst)

# Save --------------------------------------------------------------------

save(fits.long, file = outpath)
print(Sys.time() - start)