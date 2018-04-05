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
source("scripts/functions/OuterComplex.R")

omega <- 2 * pi / 24
filt.tiss <- c("WFAT")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/dat.complex.fixed_rik_genes.Robj")

f1 <- "/home/yeung/projects/nconds_results/fits_long.11_tiss_3_max.weight_raw.Robj"
# f2 <- "/home/yeung/projects/nconds_results/fits_long.11_tiss_11_max.weight_raw.Robj"
# 
# load(f1, verbose=T)
# fits.long$n.rhyth <- NULL; fits.long$n.params <- NULL
# fits.long1 <- fits.long
# load(f2, verbose=T)
# fits.long <- rbind(fits.long1, fits.long)
# rm(fits.long1)

load(f1)

# fix weight = 0 by taking lowest weight.raw
fits.best <- fits.long %>%
  group_by(gene) %>%
  filter(weight.raw == min(weight.raw))

# now to break ties of equal weight.raw, take highest weight
fits.best <- fits.best %>% 
  group_by(gene) %>%
  filter(weight == max(weight))

dat.complex <- subset(dat.complex, !tissue %in% filt.tiss)

# Summarize by top --------------------------------------------------------

fits.best <- fits.best[order(fits.best$weight, decreasing = TRUE), ]

# COLLAPSE MODEL
fits.best$model <- mapply(FilterModelByAmp, fits.best$model, fits.best$param.list, MoreArgs = list(amp.cutoff = 0.15))
fits.best$n.params <- sapply(fits.best$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.best$n.rhyth <- sapply(fits.best$model, GetNrhythFromModel)
fits.best$amp.avg <- mapply(GetAvgAmpFromParams, fits.best$param.list, fits.best$model)
fits.best$phase.sd <- mapply(GetSdPhaseFromParams, fits.best$param.list, fits.best$model)
fits.best$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.best$param.list, fits.best$model)
fits.best$phase.avg <- mapply(GePhaseFromParams, fits.best$param.list, fits.best$model)

# Get average amp and variance --------------------------------------------

head(fits.best[order(fits.best$amp.avg), ])

# Collapse models by removing models with low amps ------------------------


# high level  -------------------------------------------------------------

ggplot(fits.best, aes(x = weight)) + geom_density() + facet_wrap(~n.rhyth)
ggplot(fits.best, aes(x = weight, y = amp.avg)) + geom_point(alpha = 0.5) + facet_wrap(~n.rhyth)
ggplot(fits.best, aes(x = amp.avg)) + geom_density() + facet_wrap(~n.rhyth) + scale_x_log10()

ggplot(subset(fits.best, n.rhyth == 1), aes(x = phase.avg, y = amp.avg)) + geom_point(alpha = 0.7) + xlab("Phase (h)") + ylab("Amplitude") + facet_wrap(~model) +
  scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + 
  theme_bw(24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")

fits.tspec <- subset(fits.best, n.rhyth == 1)
# order by number of genes
fits.tspec.sum <- fits.tspec %>%
  group_by(model) %>%
  summarise(count = length(model))
fits.tspec.sum <- fits.tspec.sum[order(fits.tspec.sum$count, decreasing = TRUE), ]
fits.tspec$model <- factor(as.character(fits.tspec$model), levels = fits.tspec.sum$model)
ggplot(fits.tspec, aes(y = phase.avg, x = amp.avg)) + 
  geom_point(size=1.75) + 
  coord_polar(theta = "y") + 
  ylab("Phase (h)") + 
  xlab("Amplitude") + 
  facet_wrap(~model) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) + 
  theme_bw(24) + 
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
  expand_limits(x = 0)

m <- ggplot(data = df, aes(x = amp, y = phase, label = label)) + 
  geom_point(size = 0.5) +
  coord_polar(theta = "y") + 
  xlab(xlab) +
  ylab(ylab) +
  ggtitle(title) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))

fits.rhyth <- subset(fits.best, n.params > 0)
fits.rhyth$label <- apply(fits.rhyth, 1, function(row){
  cutoff <- 1
  if (row[8] > cutoff & row[6] > 0){  # amp.avg > cutoff only for n.rhyth > 1
    return(as.character(row[1]))  # return gene
  } else {
    return("")
  }
})
ggplot(fits.rhyth, aes(x = weight, y = amp.avg, label = label)) + geom_point(alpha = 0.25) + facet_wrap(~n.rhyth) + geom_text()
ggplot(subset(fits.rhyth, n.rhyth > 1), aes(x = phase.sd, y = amp.avg, label = label)) + geom_point(alpha = 0.25) + facet_wrap(~n.rhyth) + geom_text()
ggplot(subset(fits.rhyth, n.rhyth > 1), aes(x = phase.maxdiff, y = amp.avg, label = label)) + geom_point(alpha = 0.25) + facet_wrap(~n.rhyth) + geom_text()
ggplot(subset(fits.rhyth, n.rhyth > 1), aes(x = phase.maxdiff, y = amp.avg)) + geom_point(alpha = 0.25) + facet_wrap(~n.rhyth)


# Tissue-specific genes ---------------------------------------------------

tissues <- unique(dat.long$tissue)

for (tiss in tissues){
  fits.ts <- subset(fits.best, n.rhyth == 1 & model == tiss)
  genes.ts <- as.character(fits.ts$gene)
  print(paste("N genes for", tiss, ":", length(genes.ts)))
  if (length(genes.ts) <= length(tissues)){
    # problems with too few genes skip
    next
  }
  
  s.ts <- SvdOnComplex(subset(dat.complex, gene %in% genes.ts), value.var = "exprs.transformed")
  eigens.ts <- GetEigens(s.ts, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4)
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
  eigens.trip <- GetEigens(s.trip, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.trip$u.plot, eigens.trip$v.plot, layout = jlayout)  
}


# Aorta BFAT --------------------------------------------------------------

# load("Robjs/fits.best.collapsed_models.Robj")
# fits.best.bckup <- fits.best
fits.bfataorta <- subset(fits.best, n.rhyth == 2)
fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT|BFAT.*Aorta", fits.bfataorta$model), ]
# fits.bfataorta <- fits.bfataorta[grep("Aorta;BFAT|BFAT;Aorta", fits.bfataorta$model), ]
# fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT|BFAT.*Aorta", fits.bfataorta$model), ]

genes.bfataorta <- as.character(fits.bfataorta$gene)
print(paste("Genes in bfat-aorta:", length(genes.bfataorta)))

#outobj <- PlotHeatmapNconds(fits.bfataorta, dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5)
s.bfataorta <- SvdOnComplex(subset(dat.complex, gene %in% genes.bfataorta & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
eigens.bfataorta <- GetEigens(s.bfataorta, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.bfataorta$u.plot, eigens.bfataorta$v.plot, layout = jlayout)  


# Adr,Aorta,BFAT
fits.adrbfataorta <- subset(fits.best, n.rhyth == 3)
fits.adrbfataorta <- fits.adrbfataorta[grep("Adr.*Aorta.*BFAT", fits.adrbfataorta$model), ]

genes.adrbfataorta <- as.character(fits.adrbfataorta$gene)
print(paste("Genes in bfat-aorta:", length(genes.adrbfataorta)))

s.adrbfataorta <- SvdOnComplex(subset(dat.complex, gene %in% genes.adrbfataorta & ! tissue %in% filt.tiss), value.var = "exprs.transformed")
eigens.adrbfataorta <- GetEigens(s.adrbfataorta, period = 24, comp = 1, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens.adrbfataorta$u.plot, eigens.adrbfataorta$v.plot, layout = jlayout)  



# Pairs yo ----------------------------------------------------------------

fits.liverpairs <- subset(fits.best, n.rhyth == 3)
fits.liverpairs <- fits.liverpairs[grep("Kidney.*Liver|Liver.*Kidney", fits.liverpairs$model), ]

genes.liverpairs <- as.character(fits.liverpairs$gene)

s.liverpairs <- SvdOnComplex(subset(dat.complex, gene %in% genes.liverpairs), value.var = "exprs.transformed")
for (i in seq(3)){
  eigens.liverpairs <- GetEigens(s.liverpairs, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 2)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens.liverpairs$u.plot, eigens.liverpairs$v.plot, layout = jlayout)  
}

