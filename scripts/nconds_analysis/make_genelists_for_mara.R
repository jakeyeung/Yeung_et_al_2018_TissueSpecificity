# 2015-11-18
# Jake Yeung
# Make genelists of gene modules and run MARA.

writedir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/bic_modules"
dir.create(writedir)

# load --------------------------------------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)



# Write gene lists --------------------------------------------------------

modelname <- "Liver"

fname <- paste0(modelname, ".txt")
genes.liver <- as.character(subset(fits.best, model == modelname)$gene)

sink(file = file.path(writedir, fname))
for (g in genes.liver){
  cat(g)
  cat("\n")
}
sink()

modelname <- "TissueWide"

fname <- paste0(modelname, ".txt")
genes.tw <- as.character(subset(fits.best, n.rhyth >= 8)$gene)
sink(file = file.path(writedir, fname))
for (g in genes.tw){
  cat(g)
  cat("\n")
}
sink()


# modelname <- "AdrBFATMus"
# fname <- paste0(modelname, ".txt")
# fits.bfataorta <- fits.best[grep("Aorta.*BFAT", fits.best$model), ]
# 
# sink(file = file.path(writedir, fname))
# for (g in genes.liver){
#   cat(g)
#   cat("\n")
# }
# 
# modelname <- "AdrBFATGrep"
# fname <- paste0(modelname, ".txt")
# fits.bfataorta <- subset(fits.best, n.rhyth > 1 & n.rhyth < 11)
# fits.bfataorta <- fits.bfataorta[grep("(;|^)Aorta.*;BFAT(;|$)", fits.bfataorta$model), ]
# 
# sink(file = file.path(writedir, fname))
# for (g in genes.liver){
#   cat(g)
#   cat("\n")
# }