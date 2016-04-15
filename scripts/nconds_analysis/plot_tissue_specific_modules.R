# 2015-11-24
# plot_tissue_specific_modules.R

library(ggplot2)
library(dplyr)
library(hash)

source("scripts/functions/PlotFunctions.R")
source("scripts/functions/LongToMat.R")
# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")
load("Robjs/dat.complex.fixed_rik_genes.Robj")


# Plot SVD -------------------------------------------------------------------

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


# Plot histograms ---------------------------------------------------------


fits.tspec <- subset(fits.best, n.rhyth == 1)
# order by number of genes
fits.tspec.sum <- fits.tspec %>%
  group_by(model) %>%
  summarise(count = length(model))
fits.tspec.sum <- fits.tspec.sum[order(fits.tspec.sum$count, decreasing = TRUE), ]
fits.tspec$model <- factor(as.character(fits.tspec$model), levels = fits.tspec.sum$model)

ggplot(fits.tspec, aes(y = phase.avg, x = amp.avg)) + 
  geom_point(size=1.75, alpha = 0.25) + 
  coord_polar(theta = "y") + 
  ylab("Phase (h)") + 
  xlab("Amplitude") + 
  facet_wrap(~model) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) + 
  theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
  scale_x_continuous(limits = c(0, 3))

# plot histogram
# ggplot(subset(fits.tspec, ! model %in% c("Liver", "Adr", "BFAT", "Mus")), aes(x = phase.avg)) + 
ggplot(fits.tspec, aes(x = phase.avg)) + 
#   geom_bar(width = 1) +
  # geom_histogram(bins = 25) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~model) + 
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) +
  expand_limits(y = 0) +
  expand_limits(x = c(0, 24)) +
  xlab("Phase (h)") + 
  ylab("Count") +
  theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
  coord_polar(theta = "x")

# bin manually
# 24 and 0 are the same bin
fits.tspec$bin <- sapply(fits.tspec$phase.avg, function(p){
  if (p > 0.5){
    phase.bin <- round(p, digits = 0)
  } else {
    phase.bin <- 24
  }
  return(phase.bin)
})
fits.tspec.forhist <- fits.tspec %>%
  group_by(model, bin) %>%
  summarise(Count = length(bin))

fits.tspec.forhist$bin <- fits.tspec.forhist$bin - 0.5

# add label to indicate number of genes
n.genes.hash <- hash(as.character(fits.tspec.sum$model), fits.tspec.sum$count)
# rename tissues
fits.tspec.forhist$model.n <- as.factor(sapply(as.character(fits.tspec.forhist$model), function(m) paste0(m, " N=", n.genes.hash[[m]])))
tiss.ordered <- sapply(as.character(fits.tspec.sum$model), function(m) paste0(m, " N=", n.genes.hash[[m]]))
fits.tspec.forhist$model.n <- factor(as.character(fits.tspec.forhist$model.n), levels = tiss.ordered)

ggplot(fits.tspec.forhist, aes(x = bin, y = Count)) +
  geom_bar(stat = "identity", width = 1) + 
  facet_wrap(~model.n) + 
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) +
  expand_limits(y = 0) +
  expand_limits(x = c(0, 24)) +
  xlab("Phase (h)") + 
  ylab("Count") +
  theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
  coord_polar(theta = "x")


# plot normalized to total density
# http://wenda.io/questions/2417593/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group.html
ggplot(fits.tspec, aes(x = phase.avg)) + 
  geom_histogram(aes(y = 0.2 * ..density..), binwidth = 0.2) +  # multiply density by binwidth
  facet_wrap(~model) + 
  scale_x_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) +
  expand_limits(y = 0) +
  expand_limits(x = c(0, 24)) +
  xlab("Phase (h)") + 
  ylab("Count") +
  theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom") + 
  coord_polar(theta = "x")


# Heatmaps ----------------------------------------------------------------


PlotHeatmapNconds(subset(fits.best, model == "Liver"), dat.long, filt.tiss, jexperiment="array", blueend = -0.75, blackend = 0.75, min.n = -2.5, max.n = 2.5, remove.gene.labels = TRUE)

fits.bfataorta <- subset(fits.best, n.rhyth == 2)
fits.bfataorta <- fits.bfataorta[grep("Aorta.*BFAT|BFAT.*Aorta", fits.bfataorta$model), ]
genes.bfataorta <- as.character(fits.bfataorta$gene)
PlotHeatmapNconds(subset(fits.best, gene %in% genes.bfataorta), dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5, remove.gene.labels = TRUE)

PlotHeatmapNconds(subset(fits.best, gene %in% genes.bfataorta), subset(dat.long, tissue %in% c("Aorta", "BFAT")), filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5, jlabRow = NA, jlabCol = NA)


fits.bfataorta <- subset(fits.best, n.rhyth > 1 & n.rhyth < 11)
fits.bfataorta <- fits.bfataorta[grep("(;|^)Aorta.*;BFAT(;|$)", fits.bfataorta$model), ]
genes.bfataorta <- as.character(fits.bfataorta$gene)
PlotHeatmapNconds(subset(fits.best, gene %in% genes.bfataorta), dat.long, filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5, jlabRow = NULL, jlabCol = NA)
# PlotHeatmapNconds(subset(fits.best, gene %in% genes.bfataorta), subset(dat.long, tissue %in% c("Aorta", "BFAT")), filt.tiss, jexperiment="array", blueend = -1, blackend = 1, min.n = -2.5, max.n = 2.5, jlabRow = NA, jlabCol = NA)
