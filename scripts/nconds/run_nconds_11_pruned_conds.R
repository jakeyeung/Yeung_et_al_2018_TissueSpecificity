# Run nconds.R

setwd("/home/yeung/projects/tissue-specificity")

library("devtools")
dev_mode()
install("~/projects/ncond")  # use jake branch
library(nconds)


# Load and plot on 7 conditions -------------------------------------------

# load from this because it is easier, it is run from a code that takes a while to calculate
load("~/projects/select-rhythmic-models/data/tissue_array_dat.Robj")
load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)

conds <- c('Adr', 'Aorta', 'BFAT', 'BS', 'Cere', 'Heart', 'Hypo', 'Kidney', 'Liver', 'Lung', 'Mus')  # order matters

# scale it
library(dplyr)
library(reshape2)
dat.scaled <- subset(dat.long, experiment == "rnaseq" & tissue %in% conds) %>%
  group_by(gene, tissue) %>%
  mutate(exprs.scaled = scale(exprs, center = TRUE, scale = TRUE))

dat.filt <- dcast(dat.scaled, formula = gene ~ tissue + time, value.var = "exprs.scaled")
rownames(dat.filt) <- dat.filt$gene; dat.filt$gene <- NULL

# conds <- c('Liver', 'BFAT', 'Adr', 'Mus', 'Lung', 'Aorta', 'Kidney')
co <- length(conds)
t.vec <- as.numeric(sapply(colnames(dat.filt), function(c) strsplit(c, "_")[[1]][[2]]))

start <- Sys.time()
A <- nconds(dat.filt, out.prefix = "/home/yeung/projects/tissue-specificity/nconds_outputs/nconds_11_tissues_pruned/nconds_11_tissues", ncores = 20, max.params = 3)
print(Sys.time() - start)


dev_mode(F)

