# 2016-04-21
# Jake Yeung
# pLDA was done before with cutoff at 0.5, but no real reason to do that. Redo with 0.1!

start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")
# Load --------------------------------------------------------------------

# inf <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed_multigene/all_sites.closest.filter.genelst.dist.50000.long.bed"
# inf <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed_multigene/all_sites.closest.filter.genelst.dist.50000.long.with_peaks.bed"
inf <- "/home/yeung/data/tissue_specificity/motevo_dhs/closest_bed_multigene/all_sites.closest.filter.Mus_flat.dist.50000.long.bed"
outf <- "Robjs/N.long.Mus_flat.Robj"
N.long.livertwflat <- read.table(inf, header = FALSE, sep = "\t")  # 16 gigs
colnames(N.long.livertwflat) <- c("chromo", "start", "end", "motif", "sitecount", "gene", "dist", "peak")
save(N.long.livertwflat, file = outf)  # 2.5 gigs
print(Sys.time() - start)
