# 2015-09-29
# Jake Yeung

setwd("/home/yeung/projects/tissue-specificity")
jvalue.var <- "exprs.transformed"
jlayout <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)

outdir <- "plots/nconds/nconds_microarray_array_independent_params_11_tiss"
jtop.n <- 25

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

PlotSvdFromGeneList <- function(dat.complex, gene.list, jcomp = 1, jlabel.n = 25, jvalue.var = "exprs.transformed"){
  # Plot SVD from gene list
  s.custom <- SvdOnComplex(subset(dat.complex, gene %in% gene.list), value.var = jvalue.var)
  eigens.custom <- GetEigens(s.custom, period = 24, comp = jcomp, label.n = jlabel.n, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  multiplot(eigens.custom$u.plot, eigens.custom$v.plot, eigens.custom2$u.plot, eigens.custom2$v.plot, layout = jlayout)
  return(eigens.custom)
}

GetTopGenesFromSvd <- function(eigens.obj, top.n = 25){
  return(names(sort(eigens.obj$eigensamp, decreasing = TRUE)[1:top.n]))
}

# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/dat.complex.maxexprs4.z_nr1d1_separately.Robj")
# load("/home/yeung/projects/tissue-specificity/Robjs/nconds.rnaseq_microarray.independent.bckup.Robj", verbose = T)
# load("/home/yeung/projects/tissue-specificity/Robjs/nconds.rnaseq_microarray.independent.12_tiss.Robj", verbose = T)
load("/home/yeung/projects/tissue-specificity/Robjs/nconds.rnaseq_microarray.independent.no_WFAT.Robj", verbose = T)
filt.tiss <- c("WFAT")
# dat.long <- subset(dat.long, ! tissue %in% filt.tiss)
dat.complex <- subset(dat.complex, ! tissue %in% filt.tiss)

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


# Let's explore combos -----------------------------------------------------

fits.pairs <- subset(fits.all.max, n.tiss == 2)
fits.pairs.sum <- fits.pairs %>%
  group_by(rhyth.tiss) %>%
  summarise(count = length(rhyth.tiss))
fits.pairs.sum <- fits.pairs.sum[order(fits.pairs.sum$count, decreasing = TRUE), ]
head(fits.pairs.sum)

ggplot(OrderDecreasing(fits.pairs.sum, "rhyth.tiss", "count"), aes(x = rhyth.tiss, y = count)) + geom_bar(stat = "identity") + 
  theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1)) +
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# top pairs:
# BFAT,Liver
# Aorta,BFAT
# Adr,Liver
# Kidney,Liver
# Liver,Mus

fits.combos.sum <- subset(fits.all.max, rhyth.tiss != "") %>%  # ignore flat model
  group_by(rhyth.tiss) %>%
  summarise(count = length(rhyth.tiss)) %>%
  mutate(n.tiss = sapply(rhyth.tiss, function(s) return(length(strsplit(s, ",")[[1]]))))
fits.combos.sum <- fits.combos.sum[order(fits.combos.sum$count, decreasing = TRUE), ]
head(fits.combos.sum)


# Look at Aorta BFAT ------------------------------------------------------

# TODO: Aorta BFAT: PCA on mean to show muscle is always highly expressed: otherwise maybe there is a Aorta BFAT Mus group?

# dat.avg <- subset(dat.long, experiment == "rnaseq" & tissue %in% c("Aorta", "Heart")) %>%
#   group_by(tissue, gene) %>%
#   summarise(exprs.avg = mean(exprs))
# 
# dat.avg.diff <- dat.avg %>%
#   group_by(gene) %>%
#   summarise(abs.diff = abs(diff(exprs.avg)))
# dat.avg.diff <- dat.avg.diff[order(dat.avg.diff$abs.diff, decreasing = TRUE), ]

# fits.sub <- subset(fits.all.max, rhyth.tiss == "Aorta,BFAT")
# ggplot(fits.sub, aes(x = bicweight)) + geom_histogram(binwidth = 1/100)
# # show low bicweights
# fits.sub[order(fits.sub$bicweight), ]


# Plot combos ------------------------------------------------

# plot all combos and plot gene expression examples
tissue.combos <- as.character(fits.combos.sum$rhyth.tiss)

max.length <- 50
for (tissue.combo in tissue.combos){
  combos.genes <- subset(fits.all.max, rhyth.tiss %in% c(tissue.combo))$gene
  if (length(combos.genes) < max.length) next # too few genes are not too useful
  
  outname <- paste0(tissue.combo, ".ge", max.length, "plots.pdf")
  pdf(file.path(outdir, outname))
  print(paste(tissue.combo, "N combos:", length(combos.genes)))
  eigens <- PlotSvdFromGeneList(dat.complex, combos.genes)
  top.genes <- GetTopGenesFromSvd(eigens, top.n = jtop.n)
  for (jgene in top.genes){
    print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene)))
  }
  dev.off()
}


# Plot triplets -----------------------------------------------------------

tissue.combos <- as.character(subset(fits.combos.sum, n.tiss == 3)$rhyth.tiss)

for (tissue.combo in tissue.combos){
  combos.genes <- subset(fits.all.max, rhyth.tiss %in% c(tissue.combo))$gene
  if (length(combos.genes) < 13) next # too few genes are not too useful
  
  outname <- paste0(tissue.combo, ".triplets_plots.pdf")
  pdf(file.path(outdir, outname))
  print(paste(tissue.combo, "N combos:", length(combos.genes)))
  eigens <- PlotSvdFromGeneList(dat.complex, combos.genes)
  top.genes <- GetTopGenesFromSvd(eigens, top.n = jtop.n)
  for (jgene in top.genes){
    print(PlotGeneAcrossTissues(subset(dat.long, gene == jgene)))
  }
  dev.off()
}



# Use clustering algorithm to cluster genes with rhythmic tissues  --------

# cluster 1 to 5 using a clustering algorithm
# fits.tspec <- subset(fits.all, n.tiss >= 1 & n.tiss < tissue.wide.cutoff)
# fits.tspec$gene <- factor(as.character(fits.tspec$gene), levels = unique(fits.tspec$gene))

# cluster on all: Var.2 means flat model
bicmat.sub <- dcast(subset(fits.all, n.tiss < tissue.wide.cutoff),
                    formula = gene ~ rhyth.tiss, fill = 0, value.var = "bicweight")
rownames(bicmat.sub) <- bicmat.sub$gene; bicmat.sub$gene <- NULL

# Get clusters ------------------------------------------------------------

n.centers <- 12
clusters <- kmeans(bicmat.sub, centers = n.centers)

clusters.counts <- table(clusters$cluster)

# order by number of genes
clusters.counts <- clusters.counts[order(clusters.counts, decreasing = TRUE)]


outname <- paste0("cluster_genes_", n.centers, ".pdf")
pdf(file.path(outdir, outname))
for (i in names(clusters.counts)){
#   clusplot(bicmat.sub, clusters$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
  custom.genes <- names(clusters$cluster[which(clusters$cluster == i)])
  s.custom <- SvdOnComplex(subset(dat.complex, gene %in% custom.genes), value.var = jvalue.var)
  eigens.custom <- GetEigens(s.custom, period = 24, comp = 1, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  eigens.custom2 <- GetEigens(s.custom, period = 24, comp = 2, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  multiplot(eigens.custom$u.plot, eigens.custom$v.plot, eigens.custom2$u.plot, eigens.custom2$v.plot, layout = jlayout)
}
dev.off()
