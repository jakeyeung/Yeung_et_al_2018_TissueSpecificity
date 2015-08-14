# Jake Yeung
# 2015-08-12
# Get gene enrichment of SVD components

setwd("~/projects/tissue-specificity")

library(dplyr)
library(reshape2)

source("scripts/functions/SvdFunctions.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/ShannonEntropy.R")

RunFisher <- function(res.row){
  # From row from results of GetEnrichmentOverBins, run fisher again
  n.total <- res.row$N.genes
  x1 <- res.row$Significant  # significant and annotated
  x2 <- res.row$Annotated - x1  # not significant but annotated
  x3 <- res.row$bin - x1  # significant but not annotated
  x4 <- n.total - x3  - x2 - x1# not significant, not annotated (need total genes used to test from input)
  cont.table <- matrix(c(x1, x2, x3, x4), nrow = 2, ncol = 2, byrow = TRUE)
  fisher.out <- fisher.test(cont.table)
  # store p.value, odds ratio estimate and its 95 interval
  res.row$p.value <- fisher.out$p.value
  res.row$lower.ci <- fisher.out$conf.int[[1]]
  res.row$upper.ci <- fisher.out$conf.int[[2]]
  res.row$odds.ratio <- fisher.out$estimate
  return(res.row)
}

# Load dat ----------------------------------------------------------------

load(file = "Robjs/dat.complex.maxexprs4.Robj")
load(file = "Robjs/dat.fit.Robj")

ref.gene <- "Nr1d1"
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)

dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))

genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)

# ref gene should be max entropy as a check
dat.entropy <- dat.complex %>%
  subset(., gene %in% genes.exprs) %>%
  group_by(gene) %>%
  summarise(entropy = ShannonEntropy(Mod(exprs.adj.norm) / sum(Mod(exprs.adj.norm))))

# split genes into 3rds
N <- floor(nrow(dat.entropy) / 3)
n.genes <- nrow(dat.entropy)
low.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy), ]$gene, n = N)
high.entropy.genes <- head(dat.entropy[order(dat.entropy$entropy, decreasing = TRUE), ]$gene, n = N)
med.entropy.genes <- dat.entropy[order(dat.entropy$entropy), ]$gene[N:(n.genes - N)]

s.low.norm <- SvdOnComplex(subset(dat.complex, gene %in% low.entropy.genes), value.var = "exprs.adj")
s.high.norm <- SvdOnComplex(subset(dat.complex, gene %in% high.entropy.genes), value.var = "exprs.adj")


# Do enrichment -----------------------------------------------------------

# precalculation of these objs makes things go faster
sym2entrez <- CreateSym2Entrez()
entrez2GO <- CreateEntrez2GO()

# ilist <- c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000)
# comps <- c(1, 2, 3)
# 
# for (comp in comps){
#   gene.loadings <- names(sort(Mod(s.low.norm$u[, comp]), decreasing = TRUE))
#   fname.base <- paste0("results/GO/low_entropy/low_entropy_module", comp)
#   EnrichmentBinnedToFile(gene.loadings, ilist, fname.base, sym2entrez, entrez2GO)
# }
# 
# for (comp in comps){
#   gene.loadings <- names(sort(Mod(s.high.norm$u[, comp]), decreasing = TRUE))
#   fname.base <- paste0("results/GO/high_entropy/high_entropy_module", comp)
#   EnrichmentBinnedToFile(gene.loadings, ilist, fname.base, sym2entrez, entrez2GO)
# }


# Track circadian rhythms over bins ---------------------------------------

# ilist <- c(10, 25)

ilist <- 100 * seq(1:50)
ilist <- c(10, 25, 50, ilist)

comps <- c(1, 2, 3)
outdir <- "Robjs/GO_objs"

for (comp in comps){
  fname <- paste0("GO_high_entropy_module_", comp, ".Robj")
  outpath <- file.path(outdir, fname)
  gene.loadings <- names(sort(Mod(s.high.norm$u[, comp]), decreasing = TRUE))
  res.all <- GetEnrichmentParallel(gene.loadings, ilist, "BP", sym2entrez, entrez2GO, n.cores = 22)
  print(paste("Saving to:", outpath))
  save(res.all, file = file.path(outdir, fname))  
}

for (comp in comps){
  fname <- paste0("GO_low_entropy_module_", comp, ".Robj")
  outpath <- file.path(outdir, fname)
  gene.loadings <- names(sort(Mod(s.low.norm$u[, comp]), decreasing = TRUE))
  res.all <- GetEnrichmentParallel(gene.loadings, ilist, "BP", sym2entrez, entrez2GO, n.cores = 22)
  print(paste("Saving to:", outpath))
  save(res.all, file = file.path(outdir, fname))  
}

# AnalyzeGeneEnrichment(genes.bg = gene.loadings, genes.hit = gene.loadings[1:5000], sym2entrez, entrez2GO, FDR.cutoff = 1)
# GetEnrichment(gene.loadings, ilist[2], "BP", sym2entrez, entrez2GO)

# Find most "significant" modules -----------------------------------------

# res.all$classicFisher <- as.numeric(res.all$classicFisher)
# res.all$log10pval <- -log10(res.all$classicFisher)
# 
# res.all$bin <- factor(res.all$bin)
# 
# res.mat <- dcast(res.all, GO.ID + Term ~ bin, value.var = "log10pval")
# # res.mat <- dcast(res.all, GO.ID ~ bin, value.var = "log10pval")
# 
# rownames(res.mat) <- make.names(res.mat$Term, unique = TRUE)
# 
# res.mat <- res.mat[, 3:ncol(res.mat)]
# 
# # res.mat <- res.mat[which(rowSums(res.mat) == 0), ]
# 
# res.pca <- prcomp(t(res.mat), center = FALSE)
# 
# biplot(res.pca)

# Plot for circadian and cellular response --------------------------------
# 
# ggplot(subset(res.all, Term %in% c("response to organic substance", "single-organism metabolic process")), 
#        aes(x = bin, y = log10pval, group = Term, colour = Term, fill = Term)) + 
#   geom_point() + geom_line()

# # Run Fisher's exact test individually ------------------------------------
# 
# genes.bg <- gene.loadings
# genes.hit <- gene.loadings[1:25]
# 
# genes.bg <- ConvertSym2Entrez(genes.bg, sym2entrez)
# genes.hit <- ConvertSym2Entrez(genes.hit, sym2entrez)
# 
# # Select genes in bg that are in hit, binary matrix.
# # used in topGO new() function
# sel.genes <- factor(as.integer(genes.bg %in% genes.hit))
# names(sel.genes) <- genes.bg
# 
# # Get topGO object
# GOdata <- new("topGOdata",
#               ontology = "BP",
#               allGenes = sel.genes,
#               nodeSize = 5,
#               annot = annFUN.gene2GO,
#               gene2GO = entrez2GO)
# 
# goID <- "GO:0048511"
# gene.universe <- genes(GOdata)
# go.genes <- genesInTerm(GOdata, goID)[[1]]
# sig.genes <- sigGenes(GOdata)
# 
# my.group <- new("classicCount", testStatistic = GOFisherTest, name = "fisher",
#                 allMembers = gene.universe, groupMembers = go.genes,
#                 sigMembers = sig.genes)
# 
# runTest(my.group)
# fisher.test(contTable(my.group))
# RunFisher(res.sub[2, ], n.total = length(gene.universe))
# contTable(my.group)

