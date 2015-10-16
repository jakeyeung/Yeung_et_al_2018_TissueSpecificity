# 2015-10-12
# Jake Yeung
# Analyze nconds fit after fixing bug where I needed to fit intercepts for all RNA-Seq rather than just 1.
library(dplyr)
library(ggplot2)
library(reshape2)
library(hash)
library(gplots)

# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")


# Functions ---------------------------------------------------------------

# RangePhase <- function(phases, period = 24){
#   # get range of phase, taking modulo of period
#   
# }
# 
# AveragePhase <- function(phases, period = 24){
#   # get average phase, taking modulo 
#   
# }



# Load --------------------------------------------------------------------

load("/home/yeung/projects/nconds_results/fits_long.11_tiss_3_max.weight_raw.Robj", verbose=T)
load("Robjs/dat.long.fixed_rik_genes.Robj")

fits.raw <- fits.long

# fix weight = 0 by taking lowest weight.raw
fits.long <- fits.long %>%
  group_by(gene) %>%
  filter(weight.raw == min(weight.raw))

omega <- 2 * pi / 24

# dat.complex <- dat.long %>%
#   group_by(gene, tissue) %>%
#   do(ProjectToFrequency2(., omega, add.tissue=TRUE))
load("Robjs/dat.complex.fixed_rik_genes.Robj")

filt.tiss <- c("WFAT")

dat.complex <- subset(dat.complex, !tissue %in% filt.tiss)

fits.long$n.params <- sapply(fits.long$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long$n.rhyth <- sapply(fits.long$model, GetNrhythFromModel)

# Summarize by top --------------------------------------------------------

fits.best <- fits.long %>%
  group_by(gene) %>%
  filter(weight == max(weight))

fits.best <- fits.best[order(fits.best$weight, decreasing = TRUE), ]

# Get average amp and variance --------------------------------------------

fits.best$amp.avg <- sapply(fits.best$param.list, GetAvgAmpFromParams)

# high level  -------------------------------------------------------------

ggplot(fits.best, aes(x = weight)) + geom_density() + facet_wrap(~n.rhyth)
ggplot(fits.best, aes(x = weight, y = amp.avg, )) + geom_point(alpha = 0.1) + facet_wrap(~n.rhyth)

# Tissue-specific genes ---------------------------------------------------

fits.ts <- subset(fits.best, n.params == 1 & n.rhyth == 1)

fits.ts$param.list[[1]][["Adr.amp"]]


# PCA on tissue-specific genes --------------------------------------------

genes.ts <- as.character(fits.ts$gene)

s.ts <- SvdOnComplex(subset(dat.complex, gene %in% genes.ts), value.var = "exprs.transformed")
for (i in seq(11)){
  eigens.ts <- GetEigens(s.ts, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.ts$u.plot, eigens.ts$v.plot, layout = jlayout)  
}


# Tissue-wide genes -------------------------------------------------------

fits.tw <- subset(fits.best, n.rhyth >= 8)

genes.tw <- as.character(fits.tw$gene)

s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
for (i in seq(3)){
  eigens.tw <- GetEigens(s.tw, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)  
}


# 3param genes ------------------------------------------------------------

fits.3p <- subset(fits.best, n.params == 3)

genes.3p <- as.character(fits.3p$gene)

# write to file
sink(file = "/home/yeung/projects/nconds_results/genes_with_3_params/genes.3p.txt")
for (gene in genes.3p){
  cat(gene)
  cat("\n")
}
sink()


# Focus on pairs ----------------------------------------------------------


# Triplets ----------------------------------------------------------------


fits.trip <- subset(fits.best, n.rhyth == 3)

fits.trip.count <- fits.trip %>%
  group_by(model) %>%
  summarise(count = length(model)) %>%
  arrange(desc(count))
fits.trip.count

genes.trip <- fits.trip$gene


s.trip <- SvdOnComplex(subset(dat.complex, gene %in% genes.trip), value.var = "exprs.transformed")
for (i in seq(3)){
  eigens.trip <- GetEigens(s.trip, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.trip$u.plot, eigens.trip$v.plot, layout = jlayout)  
}


# Aorta BFAT --------------------------------------------------------------

# Aorta and BFAT
fits.bfataorta <- subset(fits.best, n.rhyth >= 2 & n.rhyth <= 4)
fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT|BFAT.*Aorta", fits.bfataorta$model), ]

genes.bfataorta <- as.character(fits.bfataorta$gene)

outobj <- PlotHeatmapNconds(fits.bfataorta, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)

s.bfataorta <- SvdOnComplex(subset(dat.complex, gene %in% genes.bfataorta & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
for (i in seq(1)){
  eigens.bfataorta <- GetEigens(s.bfataorta, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.bfataorta$u.plot, eigens.bfataorta$v.plot, layout = jlayout)  
}


# Pairs yo ----------------------------------------------------------------

fits.liverpairs <- subset(fits.best, n.rhyth == 2 | n.rhyth == 3)
fits.liverpairs <- fits.liverpairs[grep("Kidney.*Liver|Liver.*Kidney", fits.liverpairs$model), ]

genes.liverpairs <- as.character(fits.liverpairs$gene)

s.liverpairs <- SvdOnComplex(subset(dat.complex, gene %in% genes.liverpairs), value.var = "exprs.transformed")
for (i in seq(3)){
  eigens.liverpairs <- GetEigens(s.liverpairs, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.liverpairs$u.plot, eigens.liverpairs$v.plot, layout = jlayout)  
}

# Heatmap yo --------------------------------------------------------------

fits.best$model <- factor(x = fits.best$model, levels = unique(fits.best$model))
n.models <- length(unique(fits.best$model))


# Summarize by number models ----------------------------------------------

fits.best.count <- fits.best %>%
  group_by(model) %>%
  summarise(gene.count = length(gene)) %>%
  arrange(desc(gene.count))

fits.best.count.filt <- subset(fits.best.count)

fits.best.count.filt$n.params <- sapply(fits.best.count.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.best.count.filt$n.rhyth <- sapply(fits.best.count.filt$model, GetNrhythFromModel)

# Heatmap of list ---------------------------------------------------------

jmodel <- "Kidney,Liver"
# jmodel <- "BFAT,Liver"
jmodel <- "Liver"
jmodel <- "Adr"
# ref <- strsplit(jmodel, ";")[[1]][[1]]  # could be just Liver if jmodel was Liver;Kidney
# ref <- "Liver"

jmodel <- "Adr;Liver"

jmodel <- "BFAT,Liver"

jmodel <- "Kidney,Liver"
# fits.bfataorta <- subset(fits.best, n.rhyth == 2 | n.rhyth == 3)
# fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT", fits.bfataorta$model), ]
# 
# fits.best.sub <- subset(fits.best, n.rhyth <= 5)
# fits.best.sub <- fits.best.sub[grep("Liver.*Kidney|Kidney.*Liver", fits.best.sub$model), ]
# fits.best.sub <- subset(fits.best, gene %in% genes.tw)
# fits.best.sub <- subset(fits.best, n.rhyth >= 8 & n.params == 3)

jmodel <- "Aorta;BFAT"
fits.best.sub <- subset(fits.best, model == jmodel)

outobj <- PlotHeatmapNconds(fits.best.sub, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)
genes.heat <- as.character(fits.best.sub$gene)
s.heat <- SvdOnComplex(subset(dat.complex, gene %in% genes.heat & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
for (i in seq(1)){
  eigens.heat <- GetEigens(s.heat, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.heat$u.plot, eigens.heat$v.plot, layout = jlayout)  
}


# # also do some SVD to show which ones are brighter
# 
# 
# tissuenames <- sapply(colnames(outobj$mat), function(c) strsplit(c, "_")[[1]][[1]], USE.NAMES = FALSE)
# times <- sapply(colnames(outobj$mat), function(c) strsplit(c, "_")[[1]][[2]], USE.NAMES = FALSE)
# 
# mat.long <- melt(outobj$mat, value.name = "scaled.exprs", varnames = c("gene", "tissue_time"))
# mat.long$tissue <- sapply(mat.long$tissue_time, function(x) strsplit(as.character(x), "_")[[1]][[1]])
# mat.long$time <- sapply(mat.long$tissue_time, function(x) strsplit(as.character(x), "_")[[1]][[2]])
# mat.long$tissue_time <- NULL
# 
# mat.sum <- mat.long %>%
#   group_by(gene, tissue) %>%
#   summarise(abs.range = diff(range(scaled.exprs)))
# mat.sum <- mat.sum[order(mat.sum$abs.range, decreasing = TRUE), ]


# Plot clusters -----------------------------------------------------------

load("Robjs/bicmat.11_tiss_max_3.clusters.top5.bug_fixed.clusters.150.Robj", verbose=T)
load("Robjs/bicmat.11_tiss_max_3.clusters")
sort(table(clusters$cluster), decreasing = T)

clusteri <- 31
clusteri <- 16
clusteri <- 23
clusteri <- 71
clusteri <- 11
clusteri <- 15
clusteri <- 34

modelsi <- names(clusters$cluster[which(clusters$cluster == clusteri)])
fits.best.sub <- subset(fits.best, model %in% modelsi)
genes.clust <- as.character(fits.best.sub$gene)

PlotHeatmapNconds(fits.best.sub, dat.long, filt.tiss, jexperiment = "array", blueend = -0.75, blackend = 0.75, min.n = -2, max.n = 2)
s.clust <- SvdOnComplex(subset(dat.complex, gene %in% genes.clust & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
for (i in seq(1){
  eigens.clust <- GetEigens(s.clust, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.clust$u.plot, eigens.clust$v.plot, layout = jlayout)  
}


