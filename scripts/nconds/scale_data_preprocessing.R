# Jake Yeung
# 2015-10-02
# scale_data_preprocessing.R

library(dplyr)

# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")


# Plot --------------------------------------------------------------------

PlotGeneAcrossTissues(subset(dat.long, gene == "Nr1d1"))


# Scale gene by gene, tissue by tissue ------------------------------------

dat.scaled <- dat.long %>%
  group_by(gene, tissue, experiment) %>%
  mutate(exprs.scaled = scale(exprs, center = TRUE, scale = TRUE))
dat.scaled$exprs <- dat.scaled$exprs.scaled

dat.sub <- head(dat.scaled)
dat.sub$exprs <- dat.sub$exprs.scaled

# Plot --------------------------------------------------------------------

jgene <- "Rgs16"
dat.sub <- subset(dat.scaled, gene == jgene)
dat.sub$exprs <- dat.sub$exprs.scaled
PlotGeneAcrossTissues(subset(dat.sub, gene == jgene))
