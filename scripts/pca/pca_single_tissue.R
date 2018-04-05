# Do PCA on single tissue
# 2015-09-24
# Jake Yeung

# jtissue <- c("Adr", "Kidney", "Liver")
# jtissue <- c("Hypo")
# jtissue <- c("BS")
jtissue <- c("Cere")

# Define functions --------------------------------------------------------
library(wordcloud)

source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/FitRhythmic.R")

# Load data, log transform ------------------------------------------------

# dat <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt", 
#                            remove.negs = TRUE, fix.rik.xgene = TRUE)

load(file = "Robjs/dat.fit.Robj", verbose = T)
dat.fit.relamp <- GetRelampByGene(dat.fit, by.gene = ref.gene)
dat.fit.relamp <- dat.fit.relamp %>%
  group_by(gene) %>%
  mutate(max.exprs = max(int.rnaseq), min.pval = min(pval), max.relamp = max(relamp), relamp.norm = relamp / sum(relamp), amp.norm = amp / sum(amp))
genes.exprs <- unique(subset(dat.fit.relamp, max.exprs >= 4)$gene)
genes.exprs <- FixRikGenes(genes.exprs)

dat <- LoadRnaSeq()

dat <- as.matrix(dat)


# Remove non-expressed genes ----------------------------------------------

dat <- dat[which(rownames(dat) %in% genes.exprs), ]

# Optionally log2 transform -----------------------------------------------

dat <- log(dat, base = 2)

# Filter out WFAT ---------------------------------------------------------

grepstr <- paste0(jtissue, collapse = "|")
dat <- dat[, grepl(grepstr, colnames(dat))]


# Complete cases ----------------------------------------------------------

dat <- dat[complete.cases(dat), ]
dat <- dat[apply(dat, 1, function(x) !any(is.infinite(x))), ]

dat_pca <- prcomp(t(dat), center=TRUE, scale.=FALSE)

screeplot(dat_pca, type = 'l')

x <- "PC2"
y <- "PC3"
# plot tissue loadings
textplot(x = dat_pca$x[, x], y = dat_pca$x[, y], words = rownames(dat_pca$x))
# plot gene loadings
top.i <- 100
gene.loadings1 <- sort(abs(dat_pca$rotation[, x]), decreasing = TRUE)[1:top.i]
gene.loadings2 <- sort(abs(dat_pca$rotation[, y]), decreasing = TRUE)[1:top.i]
gene.loadings.names <- names(gene.loadings1)
textplot(x = gene.loadings1, y = gene.loadings2, words = gene.loadings.names)
