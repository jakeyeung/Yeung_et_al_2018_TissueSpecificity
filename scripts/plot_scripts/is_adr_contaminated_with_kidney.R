# 2015-11-27
# Jake Yeung

library(dplyr)

# Examples of contamination with kidney -----------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose=T)

dat.mean <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs))


# Find highly expressed kidney and adrenal --------------------------------

dat.mean.kidadr <- subset(dat.mean, tissue %in% c("Adr", "Kidney")) %>%
  group_by(gene) %>%
  summarise(exprs.adr = exprs.mean[1], exprs.kid = exprs.mean[2], delta = exprs.kid - exprs.adr)

# filter background genes
test <- subset(dat.mean.kidadr, exprs.adr < 4)
test <- test[order(test$delta, decreasing = T), ]

# example genes:
jgenes <- c("Umod", "Rgs16")

for (g in jgenes){
  print(PlotGeneAcrossTissues(subset(dat.long, gene == g)))
}