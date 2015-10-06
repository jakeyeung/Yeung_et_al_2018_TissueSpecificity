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



# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/NcondsAnalysisFunctions.R")

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



# Amplitude by BIC weight -------------------------------------------------

library(hash)

# do mean to trough amplitude
ampdic <- hash(paste(dat.complex$gene, dat.complex$tissue, sep = ","), 2 * Mod(dat.complex$exprs.transformed))

fits.all.max.sub <- subset(fits.all.max, gene %in% dat.complex$gene)
fits.all.max.sub$avg.amp <- apply(fits.all.max.sub, MARGIN = 1, function(row, ampdic){
  n.tiss <- row$n.tiss
  if (n.tiss == 0) return(NA)
  tissues <- strsplit(row$rhyth.tiss, ",")[[1]]  # into vector
  gene <- row$gene
  amp.vec <- rep(NA, n.tiss)
  for (i in seq(n.tiss)){
    tissue <- tissues[i]
    key <- paste(gene, tissue, sep = ",")
    amp <- ampdic[[key]]
    if (! is.null(amp)){
      amp.vec[i] <- ampdic[[key]]
    } else {
      return(NA)
    }
  }
  return(mean(amp.vec))
}, ampdic = ampdic)
head(fits.all.max.sub)

ggplot(fits.all.max.sub, aes(x = bicweight, y = avg.amp)) + geom_point(alpha = 0.30)
ggplot(fits.all.max.sub, aes(x = bicweight, y = avg.amp)) + geom_point(alpha = 0.30) + facet_wrap(~n.tiss)
ggplot(subset(fits.all.max.sub, n.tiss >= tissue.wide.cutoff), aes(x = bicweight, y = avg.amp, label = gene)) + geom_point(alpha = 0.30) + facet_wrap(~n.tiss) + geom_text()

ggplot(subset(fits.all.max.sub, n.tiss == 1), aes(x = bicweight, y = avg.amp)) + geom_point(alpha = 0.3)
ggplot(subset(fits.all.max.sub, n.tiss == 2), aes(x = bicweight, y = avg.amp)) + geom_point(alpha = 0.3)



# Look at sum of bicweights -----------------------------------------------

fits.all.sum <- fits.all %>%
  group_by(gene) %>%
  summarise(bicweight.mean = mean(bicweight), bicweight.var = var(bicweight), bicweight.max = max(bicweight), bicweight.sum = sum(bicweight))

fits.max.dic <- hash(fits.all.max$gene, fits.all.max$rhyth.tiss)
fits.all.sum$rhyth.tiss <- sapply(as.character(fits.all.sum$gene), function(jgene){
  rhyth.tiss <- fits.max.dic[[jgene]]
  if (is.null(rhyth.tiss)){
    return(NA)
  } else {
    return(rhyth.tiss)
  }
})

ggplot(fits.all.sum, aes(x = bicweight.sum, bicweight.var)) + geom_point(alpha = 0.25)

fits.all.sum[order(fits.all.sum$bicweight.var), ]


# What are distribution of flat set ---------------------------------------

ggplot(subset(fits.all.max, n.tiss == 6), aes(x = bicweight)) + geom_density()
ggplot(fits.all.max, aes(x = bicweight)) + geom_histogram(binwidth = 1/100) + facet_wrap(~n.tiss)
ggplot(subset(fits.all.max, n.tiss > 0), aes(x = bicweight)) + geom_histogram(binwidth = 1/100) + facet_wrap(~n.tiss)

# Look at individual combos -----------------------------------------------

# do PCA on means
dat.means <- subset(dat.long, experiment == "rnaseq") %>%
  group_by(gene, tissue, experiment) %>%
  summarise(exprs.mean = mean(exprs))


# Aorta-BFAT --------------------------------------------------------------


tissue.comb <- "Aorta,BFAT"

fits.comb <- subset(fits.all.max, rhyth.tiss == tissue.comb)
genes.comb <- fits.comb$gene
PlotTissuePcaFromGenelist(dat.means, genes.comb)

# BFAT-Mus ----------------------------------------------------------------

tissue.comb <- "BFAT,Mus"

fits.comb <- subset(fits.all.max, rhyth.tiss == tissue.comb)
genes.comb <- fits.comb$gene
PlotTissuePcaFromGenelist(dat.means, genes.comb)

# BFAT-Aorta-Mus ----------------------------------------------------------

tissue.comb <- "Aorta,BFAT,Mus"

fits.comb <- subset(fits.all.max, rhyth.tiss == tissue.comb)
genes.comb <- fits.comb$gene
PlotTissuePcaFromGenelist(dat.means, genes.comb)

# Liver-Kidney ------------------------------------------------------------

tissue.comb <- "Kidney,Liver"

fits.comb <- subset(fits.all.max, rhyth.tiss == tissue.comb)
genes.comb <- fits.comb$gene
PlotTissuePcaFromGenelist(dat.means, genes.comb)

# Muscle ------------------------------------------------------------------

tissue.comb <- "Mus"

fits.comb <- subset(fits.all.max, rhyth.tiss == tissue.comb)
genes.comb <- fits.comb$gene
PlotTissuePcaFromGenelist(dat.means, genes.comb)
eigens <- PlotSvdFromGeneList(dat.complex, genes.comb)


# Liver -------------------------------------------------------------------

tissue.comb <- "Liver"

fits.comb <- subset(fits.all.max, rhyth.tiss == tissue.comb)
genes.comb <- fits.comb$gene
PlotTissuePcaFromGenelist(dat.means, genes.comb)
eigens <- PlotSvdFromGeneList(dat.complex, genes.comb)

# Cluster by BIC weights --------------------------------------------------

# cluster on all: Var.2 means flat model


bicmat.sub <- dcast(subset(fits.all), fun.aggregate = max,
                    formula = gene ~ rhyth.tiss, fill = 0, value.var = "bicweight")
rownames(bicmat.sub) <- bicmat.sub$gene; bicmat.sub$gene <- NULL
dim(bicmat.sub)
bicmat.sub <- bicmat.sub[complete.cases(bicmat.sub), ]
dim(bicmat.sub)

# Get clusters ------------------------------------------------------------

load("Robjs/bicmat.sub.Robj", verbose = T)
# d <- dist(t(bicmat.sub), method = "euclidean") # distance matrix
d <- dist(bicmat.sub, method = "euclidean") # distance matrix
clustfit <- hclust(d, method="ward") 
plot(clustfit) # display dendogram
groups <- cutree(clustfit, k=40) # cut tree into 120 clusters
# draw dendogram with red borders around the 5 clusters 
# rect.hclust(clustfit, k=50, border="red")
# save(d, clustfit, file = "Robjs/nconds.cluster_2048_models.Robj")
# save(d, clustfit, file = "Robjs/nconds.cluster_2048_models.by_gene.Robj")


# Plot genes in clustered models (by gene) --------------------------------

load(file = "Robjs/nconds.cluster_2048_models.by_gene.Robj", verbose = T)
groups <- cutree(clustfit, k = 450)
gene.cluster <- "Trim2"
clusteri <- groups[which(names(groups) == gene.cluster)]
cluster.name <- names(groups[which(groups == clusteri)])
cluster.genes <- subset(fits.all.max, gene %in% cluster.name)$gene
eigens <- PlotSvdFromGeneList(dat.complex, cluster.genes)


# Look at some genes ------------------------------------------------------

jgene <- "Trim2"
subset(fits.all, gene == jgene)
PlotGeneAcrossTissues(subset(dat.long, gene == jgene))

# Plot genes in clustered models ------------------------------------------

# find tissue-wide model
# Adr,Aorta,BFAT,BS,Cere,Heart,Hypo,Kidney,Liver,Lung,Mus
# tissue.combo <- "Adr,Aorta,BFAT,BS,Cere,Heart,Hypo,Kidney,Liver,Lung,Mus"
# tissue.combo <- "Liver"
clusteri <- groups[which(names(groups) == tissue.combo)]

cluster.name <- names(groups[which(groups == clusteri)])

# Plot genes in these clusters --------------------------------------------

cluster.genes <- subset(fits.all.max, rhyth.tiss %in% cluster.name)$gene
eigens <- PlotSvdFromGeneList(dat.complex, cluster.genes)



# n.centers <- 250
# clusters <- kmeans(bicmat.sub, centers = n.centers, iter.max = 100)
# 
# clusters.counts <- table(clusters$cluster)
# 
# 
# # order by number of genes
# clusters.counts <- clusters.counts[order(clusters$withinss, decreasing = TRUE)]
# 
# 
# 
# # Observe clusters --------------------------------------------------------
# 
# clusters.counts
# clusters$withinss
# 
# clusteri <- 168
# custom.genes <- names(clusters$cluster[which(clusters$cluster == clusteri)])
# s.custom <- SvdOnComplex(subset(dat.complex, gene %in% custom.genes), value.var = jvalue.var)
# eigens.custom <- GetEigens(s.custom, period = 24, comp = 1, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
# eigens.custom2 <- GetEigens(s.custom, period = 24, comp = 2, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
# multiplot(eigens.custom$u.plot, eigens.custom$v.plot, eigens.custom2$u.plot, eigens.custom2$v.plot, layout = jlayout)
