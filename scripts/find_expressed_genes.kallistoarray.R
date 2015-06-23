# 2015-06-22
# Jake Yeung
# Find genes expressed above a threshold

# constants for handlign Log2 transformation
scale.factor <- 100
pseudocount <- 1

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadLong.R")
source("scripts/functions/MixtureModelFunctions.R")
source("scripts/functions/WriteGeneListMat.R")

# Load --------------------------------------------------------------------

ka.long <- LoadLong(scale.factor = 100, pseudocount = 1)
kallisto.wide <- LoadKallistoGene(rna.seq.path, form = "wide")
kallisto.wide <- log2(scale.factor * kallisto.wide + pseudocount)

exprs.log2 <- unlist(as.matrix(kallisto.wide))
exprs.log2 <- exprs.log2[which(exprs.log2 > 0)]
# exprs.log2 <- log2(scale.factor * exprs.norm + pseudocount)

# Plot density  -----------------------------------------------------------

plot(density(exprs.log2))


# Fit mixture model -------------------------------------------------------

cutoff <- FindCutoff(x = exprs.log2, lambdas = c(0.2, 0.8), mus = c(1, 10), k = 2)  # max=5.572458


# Filter genes based on median gene expression ----------------------------

row.meds <- apply(kallisto.wide, 1, function(x) median(x))

expressed.genes <- names(row.meds[which(row.meds >= cutoff)])


# Filter further for merged array and rnaseq ------------------------------

ka.sub <- subset(ka.long, gene %in% expressed.genes)

ka.sub <- ka.sub %>%
  group_by(gene) %>%
  mutate(exprs.centered = scale(exprs, center = TRUE, scale = FALSE))

tissues <- unique(as.character(ka.sub$tissue))
outdir="/home/yeung/projects/tissue-specificity/data/gene_lists/expressed_genes_kallisto"
dir.create(outdir)

WriteGeneListMat2(tissues = tissues, genelist = expressed.genes, dat = ka.sub, valvar = "exprs.centered", outdir = outdir)

