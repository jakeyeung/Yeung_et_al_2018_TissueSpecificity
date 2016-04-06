# 2015-04-05
# Jake Yeung
# read_filter_motif_files.Rj
# do for multigene list

library(data.table)
start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/SitecountsFunctions.R")

library(parallel)

filtered_closest_dir <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bd_multigene_by_motif_filter"

jfiles <- list.files(filtered_closest_dir, full.names = TRUE)
N.long.filt <- mclapply(jfiles, function(f){
  return(ReadSitecountsMotif(path = f, cnames = c("chromo", "start", "end", "motif_peak", "sitecount", "gene", "dist")))
}, mc.cores = 19)
# print(N.long.filt)
N.long.filt <- rbindlist(N.long.filt)
# N.long.filt <- do.call(rbind, N.long.filt)

# save(N.long.filt, file = "Robjs/N.long.filt.liver_genes.all_motifs.Robj")
# N.long.filt <- subset(N.long.filt)
save(N.long.filt, file = "Robjs/N.long.multigene.distfilt.50000.Robj")

print(Sys.time() - start)
