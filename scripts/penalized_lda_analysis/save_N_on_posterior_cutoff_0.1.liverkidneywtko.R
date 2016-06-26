# 2015-04-05
# Jake Yeung
# read_filter_motif_files.Rj
# do for multigene list

start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")
# Load --------------------------------------------------------------------

inf <- "/home/yeung/data/tissue_specificity/motevo_dhs/liver_kidney_wtko_modules/all_sites.closest.filter.genelst.dist.50000.long.bed"
outf <- "Robjs/liver_kidney_atger_nestle/N.long.3wtmodules.Robj"
N.long.livertwflat <- read.table(inf, header = FALSE, sep = "\t")  # 16 gigs
colnames(N.long.livertwflat) <- c("chromo", "start", "end", "motif", "sitecount", "gene", "dist", "peak")
save(N.long.livertwflat, file = outf)  # 2.5 gigs
print(Sys.time() - start)