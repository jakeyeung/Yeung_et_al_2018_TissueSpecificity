# 2015-09-29
# Jake Yeung

library(dplyr)
library(ggplot2)
library(reshape2)
library(cluster) 

# Functions ---------------------------------------------------------------

GetTissuesFromCoefFit <- function(fit.coef.names, tissue.colname = "tissue"){
  # Get tissues involved by grepping everything after tissue.colname
  match.names.i <- grepl(tissue.colname, fit.coef.names) & !grepl(":", fit.coef.names)
  match.names <- fit.coef.names[match.names.i]
  tissues <- sapply(match.names, function(jname) strsplit(jname, split = tissue.colname)[[1]][[2]], USE.NAMES = FALSE)
  return(tissues)
}

ExtraParamsFromFit <- function(fit.coef, tissues){
  # From coefficients of fit, extract parameters
  # colnames should contain tissue, RNA-Seq intercept, phase, amplitude, and BIC
  cnames <- tissues
  # TODO
}

SubsetByMaxBicWeight <- function(dat){
  max.i <- which.max(dat$bicweight)
  return(dat[max.i, ])
}

# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/dat.complex.maxexprs4.z_nr1d1_separately.Robj")
# load("/home/yeung/projects/tissue-specificity/Robjs/nconds.rnaseq_microarray.independent.bckup.Robj", verbose = T)
load("/home/yeung/projects/tissue-specificity/Robjs/nconds.rnaseq_microarray.independent.12_tiss.Robj", verbose = T)

tissues <- GetTissuesFromCoefFit(names(fits.all$mod[[1]]$fit))

rhyth.tiss.bckup <- fits.all$rhyth.tiss
fits.all$rhyth.tiss <- sapply(rhyth.tiss.bckup, function(x) return(x[[1]]))
fits.all$n.tiss <- sapply(fits.all$rhyth.tiss, function(s) return(length(strsplit(s, ",")[[1]])))


# Filter by max -----------------------------------------------------------

fits.all.max <- fits.all %>%
  group_by(gene) %>%
  do(SubsetByMaxBicWeight(.))


# Plot distribution of models
fits.counts <- fits.all.max %>%
  group_by(rhyth.tiss) %>%
  summarise(count = length(rhyth.tiss)) %>%
  mutate(n.tiss = sapply(rhyth.tiss, function(s) return(length(strsplit(s, ",")[[1]]))))

fits.counts$rhyth.tiss <- as.factor(fits.counts$rhyth.tiss)
fits.counts <- OrderDecreasing(fits.counts, "rhyth.tiss", "count")
fits.counts <- fits.counts[order(fits.counts$count, decreasing = TRUE), ]


# Cluster all genes greater than or equal to 6 as one cluster -------------

tissue.wide.cutoff <- 6
tissue.wide.genes <- subset(fits.all.max, n.tiss >= 6)$gene
print(paste("Tissue wide gene count:", length(tissue.wide.genes)))
s.custom <- SvdOnComplex(subset(dat.complex, gene %in% tissue.wide.genes), value.var = jvalue.var)
eigens.custom <- GetEigens(s.custom, period = 24, comp = 1, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
eigens.custom2 <- GetEigens(s.custom, period = 24, comp = 2, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
multiplot(eigens.custom$u.plot, eigens.custom$v.plot, eigens.custom2$u.plot, eigens.custom2$v.plot, layout = jlayout)


# Use clustering algorithm to cluster genes with rhythmic tissues  --------

# cluster 1 to 5 using a clustering algorithm
# fits.tspec <- subset(fits.all, n.tiss >= 1 & n.tiss < tissue.wide.cutoff)
# fits.tspec$gene <- factor(as.character(fits.tspec$gene), levels = unique(fits.tspec$gene))

# cluster on all: Var.2 means flat model
bicmat.sub <- dcast(subset(fits.all, n.tiss < tissue.wide.cutoff),
                    formula = gene ~ rhyth.tiss, fill = 0, value.var = "bicweight")
rownames(bicmat.sub) <- bicmat.sub$gene; bicmat.sub$gene <- NULL

# Get clusters ------------------------------------------------------------

n.centers <- 25
clusters <- kmeans(bicmat.sub, centers = n.centers)

clusters.counts <- table(clusters$cluster)

# order by number of genes
clusters.counts <- clusters.counts[order(clusters.counts, decreasing = TRUE)]

jvalue.var <- "exprs.transformed"
jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
outdir <- "plots/nconds/nconds_microarray_array_independent_params_12_tiss"
outname <- paste0("cluster_genes_", n.centers, ".pdf")
pdf(file.path(outdir, outname))
for (i in names(clusters.counts)){
  clusplot(bicmat.sub, clusters$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
  custom.genes <- names(clusters$cluster[which(clusters$cluster == i)])
  s.custom <- SvdOnComplex(subset(dat.complex, gene %in% custom.genes), value.var = jvalue.var)
  eigens.custom <- GetEigens(s.custom, period = 24, comp = 1, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  eigens.custom2 <- GetEigens(s.custom, period = 24, comp = 2, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  multiplot(eigens.custom$u.plot, eigens.custom$v.plot, eigens.custom2$u.plot, eigens.custom2$v.plot, layout = jlayout)
}
dev.off()



# pdf("/home/yeung/projects/tissue-specificity/plots/nconds/nconds_microarray_rnaseq_12_tissues_heatmap/heatmap.pdf")
# PlotRelampHeatmap(bicmat.sub, jtitle = "heatmap", blackend = 0.5, yellowstart = 0.51, maxval = 1, dist.method = "euclidean")
# dev.off()

# ggplot(subset(fits.counts, n.tiss > 1), aes(x = rhyth.tiss, y = count)) + 
#   geom_bar(stat = "identity") + xlab("") + theme_bw(24) +
#   theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) +
#   scale_x_discrete(labels = "")
# 
# bicmat <- dcast(fits.all, formula = gene ~ rhyth.tiss, fun.aggregate = max, fill = 0, value.var = "bicweight")
# rownames(bicmat) <- bicmat$gene; bicmat$gene <- NULL
# 
# # remove rows where Var.2 (flat model) is > 0
# bicmat <- bicmat[which(bicmat$Var.2 == 0), ]; bicmat$Var.2 <- NULL
# 
# PlotRelampHeatmap(bicmat, jtitle = "heatmap", blackend = 0.5, yellowstart = 0.51, maxval = 1, dist.method = "manhattan")

