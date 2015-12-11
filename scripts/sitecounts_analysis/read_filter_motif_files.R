# 2015-12-09
# Jake Yeung
# read_filter_motif_files.R
# read each motif file and filter it for gene list, then concatenate into a larger file and save to output

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/SitecountsFunctions.R")

library(parallel)

# Load --------------------------------------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)

# to change
jmodel <- ""
outfile <- "Robjs/N.long.flat_genes.all_motifs.Robj"

liver.genes <- subset(fits.best, model == jmodel)$gene

print(liver.genes)

sitecounts.dir <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed"

jfiles <- list.files(sitecounts.dir)

N.long <- mclapply(jfiles, function(f){
  fpath <- file.path(sitecounts.dir, f)
  return(subset(ReadSitecountsMotif(path = fpath), gene %in% liver.genes))
}, mc.cores = 19)

N.long <- do.call(rbind, N.long)
# save(N.long, file = "Robjs/N.long.liver_genes.all_motifs.Robj")
save(N.long, file = outfile)
