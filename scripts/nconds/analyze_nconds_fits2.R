# 2015-10-12
# Jake Yeung
# Analyze nconds fit after fixing bug where I needed to fit intercepts for all RNA-Seq rather than just 1.
library(dplyr)
library(ggplot2)
library(reshape2)
library(hash)
# Source ------------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")

# Load --------------------------------------------------------------------

load("/home/yeung/projects/nconds_results/fits_long.11_tiss_3_max.weight_raw.Robj", verbose=T)
load("Robjs/dat.long.fixed_rik_genes.Robj")

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
ggplot(fits.best, aes(x = weight, y = amp.avg)) + geom_point(alpha = 0.1) + facet_wrap(~n.rhyth)

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

fits.tw <- subset(fits.best, n.params == 3 & n.rhyth >= 8)

genes.tw <- as.character(fits.tw$gene)

s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
for (i in seq(3)){
  eigens.tw <- GetEigens(s.tw, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)  
}


# 3param genes ------------------------------------------------------------

fits.3p <- subset(fits.best, n.params == 3 & weight > 0)

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
fits.bfataorta <- subset(fits.best, n.rhyth == 2 | n.rhyth == 3)
fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT", fits.bfataorta$model), ]

genes.bfataorta <- as.character(fits.bfataorta$gene)

s.bfataorta <- SvdOnComplex(subset(dat.complex, gene %in% genes.bfataorta), value.var = "exprs.transformed")
for (i in seq(3)){
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

fits.best.count.filt <- subset(fits.best.count, gene.count > 1)


# Heatmap of list ---------------------------------------------------------

jmodel <- "Kidney,Liver"
jmodel <- "Kidney;Liver"
ref <- strsplit(jmodel, ";")[[1]][[1]]  # could be just Liver if jmodel was Liver;Kidney
fits.best.sub <- subset(fits.best, model == jmodel)
genes <- as.character(fits.best.sub$gene)

dat.sub <- subset(dat.long, gene %in% genes & experiment == "array" & !tissue %in% filt.tiss)

# center and scale
dat.sub <- dat.sub %>%
  group_by(gene, tissue) %>%
  mutate(exprs.scaled = scale(exprs, center = TRUE, scale = TRUE))


mat <- dcast(dat.sub, gene ~ tissue + time, value.var = "exprs.scaled")
rownames(mat) <- mat$gene; mat$gene <- NULL
head(mat)

# sort by phases
phases.dic.keys <- as.character(fits.best.sub$gene)
phases.dic.vals <- sapply(fits.best.sub$param.list, function(p) p[[paste0(ref, ".phase")]])
phases.dat <- data.frame(gene = phases.dic.keys, phase = phases.dic.vals)
phases.dat <- phases.dat[order(phases.dat$phase), ]

# sort by phases
mat <- mat[as.character(phases.dat$gene), ]
head(mat)

# heatmap
n.co <- length(unique(dat.sub$tissue))
time <- rep(unique(dat.sub$time), n.co)
condi_name <- as.character(unique(dat.sub$tissue))
cond_begins <- seq(1, length(time), length(time) / n.co)# represent index start of next condition. 
cond_mid <- mean(c(1, length(time) / n.co))  # mid distance between first sample in cond1 and last sample in cond1.
VLINE_X_LOC <<- cond_begins - 0.5  # Subtract 0.5 to get it between two samples.
TEXT_X_LOC <<- cond_begins + cond_mid  # a vector in middle of sample, can put text labels conveniently.
TEXT_Y_LOC <<- 0.99 * length(genes)  # number of genes
C_NAME <<- unique(condi_name)

heatmap(as.matrix(mat), 
        Rowv = NA, 
        Colv = NA, 
        ylab = NA ,
        labCol = paste('CT',time, sep = "_"),
        labRow = NA,
        scale = NULL,
        add.expr = c(abline(v = VLINE_X_LOC, 
                            col = 'white'),
                     text(x=TEXT_X_LOC, 
                          y=TEXT_Y_LOC, 
                          labels = C_NAME,
                          col = 'white')),
        main = paste("model #Genes",sep =" "), 
        col = colorRampPalette(c('blue','black','red'))(1000))

