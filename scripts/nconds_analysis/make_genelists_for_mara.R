# 2015-11-18
# Jake Yeung
# Make genelists of gene modules and run MARA.

writedir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/bic_modules"
dir.create(writedir)

# load --------------------------------------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.20.phase_sd_maxdiff_avg.Robj", verbose=T)



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

modelname <- "AdrBFATMus"
fname <- paste0(modelname, ".txt")
fits.bfataorta <- fits.best[grep("Aorta.*BFAT", fits.best$model), ]

sink(file = file.path(writedir, fname))
for (g in genes.liver){
  cat(g)
  cat("\n")
}