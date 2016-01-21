# 2015-12-09
# Jake Yeung
# read_filter_motif_files.R
# read each motif file and filter it for gene list, then concatenate into a larger file and save to output

library(data.table)
start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/SitecountsFunctions.R")

library(parallel)

filtered_closest_dir <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed_filtered"

jfiles <- list.files(filtered_closest_dir, full.names = TRUE)
N.long.filt <- mclapply(jfiles, function(f){
  return(ReadSitecountsMotif(path = f))
}, mc.cores = 19)
# print(N.long.filt)
N.long.filt <- rbindlist(N.long.filt)
# N.long.filt <- do.call(rbind, N.long.filt)

# save(N.long.filt, file = "Robjs/N.long.filt.liver_genes.all_motifs.Robj")
N.long.filt <- subset(N.long.filt, dist <= 1000)
save(N.long.filt, file = "Robjs/N.long.all_genes.all_signif_motifs.Robj")

print(Sys.time() - start)
