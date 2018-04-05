# 2015-10-05
# Jake Yeung
# save_chunks_run_nconds.R

library(dplyr)

setwd("/home/yeung/projects/tissue-specificity")


# Source ------------------------------------------------------------------

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")

dat.long.path <- "Robjs/dat.long.fixed_rik_genes.Robj"


# Load --------------------------------------------------------------------

load(dat.long.path)

dat.long <- subset(dat.long, tissue != "WFAT")

tissues <- as.character(unique(dat.long$tissue))

# For expressed genes only ------------------------------------------------

# for expressed genes only
dat.exprs <- dat.long %>%
  filter(experiment == "rnaseq") %>%
  group_by(tissue, gene) %>%
  summarise(exprs.mean = mean(exprs)) %>%
  filter(exprs.mean > 5)

genes.exprs <- unique(dat.exprs$gene)

dat.sub <- subset(dat.long, tissue %in% tiss.test & gene %in% genes.exprs.sub)
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = tiss.test)
dat.sub$gene <- factor(as.character(dat.sub$gene), levels = genes.exprs.sub)

ChunkDatGenesToFile(dat.sub, write.dir = "/home/yeung/projects/tissue-specificity/data/datlong_chunks_by_gene/expressed_genes")
