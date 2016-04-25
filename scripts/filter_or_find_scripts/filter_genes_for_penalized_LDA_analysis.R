# 2016-04-20
# Jake Yeung
# filter_genes_for_penalized_LDA_analysis.R
# Filter genes to limit size of motif bed files without compromising with cutoffs
# use this to run 


# load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)


# Get Liver, flat, tissue-wide --------------------------------------------

liver.genes <- as.character(subset(fits.best, model == "Liver")$gene)
tw.genes <- as.character(subset(fits.best, n.rhyth >= 8)$gene)
flat.genes <- as.character(subset(fits.best, model == "")$gene)


# Write to output ---------------------------------------------------------

genes.list <- c(liver.genes, tw.genes, flat.genes)

sink(file = "data/gene_lists/liverrhyth_tissuewide_flatgenes.txt")
for (g in genes.list){
  cat(g); cat("\n")
}
sink()