# 2015-09-07
# Jake Yeung
# find_expressed_genes_by_tissue.deseq.R

library(ggplot2)
library(mixtools)
library(dplyr)

outdir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/expressed_genes_deseq_by_tissue"

# Sources -----------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/MixtureModelFunctions.R")


# Load --------------------------------------------------------------------


dat.long <- LoadArrayRnaSeq(fix.rik.xgene = TRUE)


# Show plot ---------------------------------------------------------------


ggplot(subset(dat.long, experiment == "rnaseq" & exprs > 0), aes(x = exprs)) + geom_density() + facet_wrap(~tissue)



# Fit ---------------------------------------------------------------------


exprs.vec <- subset(dat.long, experiment == "rnaseq" & exprs > 0)$exprs

cutoff <- FindCutoff(exprs.vec, lambdas = c(0.2, 0.8), mus = c(1, 10))


# Filter out each gene for each tissue ------------------------------------

tissues <- unique(as.character(dat.long$tissue))

dat.max <- dat.long %>%
  subset(., experiment == "rnaseq") %>%
  group_by(tissue, gene) %>%
  summarise(exprs.max = max(exprs))
head(dat.max)

for (tiss in tissues){
  # get list of genes above cutoff
  dat.sub <- subset(dat.max, tissue == tiss & exprs.max >= cutoff$maximum)
  # write to file
  sink(file = file.path(outdir, paste0(tiss, ".genelist.txt")))
  for (g in unique(dat.sub$gene)){
    cat(g)
    cat("\n")
  }
  sink()
}

# 
# # fit gamma
# gscale <- 20
# gmeans <- c(1, 11)
# glambdas <- c(0.2, 0.8)
# mixmdl <- gammamixEM(exprs.vec, lambda = glambdas, alpha = c(gscale, gscale), beta = c(gmeans[1] / gscale, gmeans[2] / gscale))
# 
# cutoff.gamma.log2 <- PlotGammaMixmdl(mixmdl, exprs.vec)
