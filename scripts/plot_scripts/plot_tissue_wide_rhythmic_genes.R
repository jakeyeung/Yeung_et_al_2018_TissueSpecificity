# 2016-01-22
# Jake Yeung
# Plot tissue wide rhythmic genes

setwd("/home/yeung/projects/tissue-specificity")
source("scripts/functions/PlotGeneAcrossTissues.R")

# Load --------------------------------------------------------------------

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)
fits.sub <- subset(fits.best, n.rhyth >= 8 & amp.avg > 0.3)

(gene.lst <- fits.sub$gene[order(fits.sub$amp.avg, decreasing = TRUE)])

dat.sub <- subset(dat.long, gene %in% gene.lst)

# pdf("plots/gene_expressions/tissue_wide_rhythmic_genes.pdf")
# for (g in gene.lst){
#   print(PlotGeneAcrossTissues(subset(dat.sub, gene == g)))
# }
# dev.off()

sink(file = "textfiles/tissue_wide_rhythmic_genes.txt")
for (g in gene.lst){
  cat(g)
  cat("\n")
}
sink()