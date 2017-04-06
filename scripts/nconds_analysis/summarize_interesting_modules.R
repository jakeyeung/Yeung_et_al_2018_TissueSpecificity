# 2017-01-13
# Jake Yeung
# summarize_interesting_modules.R

rm(list=ls())

library(dplyr)
library(reshape2)
library(gplots)

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
# source("/home/yeung/projects/tissue-specificity/scripts/functions/PhaseColorFunctions.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/ColorFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/SvdFunctions.R")


# Load --------------------------------------------------------------------

# load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid
load("Robjs/dat.complex.fixed_rik_genes.Robj", v=T); dat.complex <- subset(dat.complex, tissue != "WFAT")

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
dat.wtko$tissue <- paste(dat.wtko$tissue, as.character(dat.wtko$geno), sep = "")
dat.wtko$tissue <- gsub("SV129", "", dat.wtko$tissue)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch
load("Robjs/N.long.promoters_500.Robj", v=T)

load("Robjs/fits.relamp.Robj", v=T)

# Summarize top hogenesch -------------------------------------------------

subset(fits.best, n.rhyth > 1 & n.rhyth < 8) %>% group_by(model) %>% summarise(n.models = length(gene)) %>% arrange(desc(n.models))


fits.best.sub <- subset(fits.best, n.rhyth >= 8)
fits.best.sub <- subset(fits.best, model == "Kidney,Liver")
fits.best.sub <- subset(fits.best, model == "Liver")
fits.best.sub <- subset(fits.best, model %in% c("Kidney;Liver", "Kidney,Liver"))

fits.lst <- list(subset(fits.best, model == "Liver"), 
                 subset(fits.best, model %in% c("BFAT;Mus")), 
                 subset(fits.best, model %in% c("BFAT,Mus")), 
                 subset(fits.best, model %in% c("Heart,Liver")), 
                 subset(fits.best, model %in% c("Kidney,Liver")), 
                 subset(fits.best, model %in% c("Kidney;Liver")), 
                 subset(fits.best, model %in% c("Aorta,BFAT,Liver")),
                 subset(fits.best, model %in% c("Liver;Aorta,BFAT")),
                 subset(fits.best, n.rhyth >= 8))
jgenes <- as.character(unique(fits.best.sub$gene))
dat.sub <- subset(dat.long, gene %in% jgenes)
# filt.tiss <- c("WFAT", "Cere", "Hypo", "BS", "Mus", "Aorta", "Heart", "Adr", "BFAT", "Lung")
filt.tiss <- c("WFAT")


jblueend <- -0.5
jblackend <- 0.5

jblueend <- -0.25
jblackend <- 0.25
# jblueend <- -0.2
# jblackend <- 0.2

jmin.n <- -1.5
jmax.n <- 1.5

# reorder dat.long tissues
dat.long$tissue <- factor(as.character(dat.long$tissue), levels = c("Liver", "Kidney", "BFAT", "Lung", "Adr", "Aorta", "Heart", "Mus", "Cere", "BS", "Hypo", "WFAT"))


# Check HIC1 sitecounts ---------------------------------------------------

pdf(paste0("plots/heatmaps/polarplots_tissuewide_hic1_sitecounts.pdf"))
N.hic <- subset(N.long, motif == "HIC1.p2")
size.hash <- hash(as.character(N.hic$gene), 2 * N.hic$sitecount)

fits <- fits.lst[[9]]  # tissuewide

# color if greater than some threshold
jbins <- 300
jcolors <- colorRampPalette(brewer.pal(9, "Greys"))(jbins)
site.min <- 0
site.max <- 1.5
x <- seq(from = site.min, to = site.max, length.out = jbins)
# assign sitecount to hash by assigning to nearest value
jcolors.i <- sapply(N.hic$sitecount, function(s) which.min(abs(s - x)))
col.hash <- hash(as.character(N.hic$gene), jcolors[jcolors.i])
# size.hash <- hash(as.character(N.hic$gene), 2*x[jcolors.i])
# col.hash <- hash(as.character(N.hic$gene), ifelse(N.hic$sitecount > 1.5, "blue", "gray"))

s <- SvdOnComplex(subset(dat.complex, gene %in% genes), value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = comp, label.n = Inf, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, jsize = 16, dotsize=size.hash, dot.col = col.hash, disable.repel = TRUE, disable.text = TRUE)
eigens.genenames <- GetEigens(s, period = 24, comp = comp, label.n = 50, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, jsize = 16, dotsize=size.hash , disable.repel = FALSE, disable.text = FALSE)
print(eigens$u.plot + ylab("ZT") + ggtitle(paste0("MinAndMax:", paste(signif(range(N.hic$sitecount), digits = 2), collapse = ","))))
dev.off()

# Plot modules ------------------------------------------------------------


  
pdf(paste0("plots/heatmaps/polarplots_heatmaps_hogenesch_model_selection.blacklims", jblueend, jblackend, ".minmax.", jmin.n, jmax.n, ".pdf"))
lapply(fits.lst, function(fits){
  print(PlotHeatmapNconds(fits, dat.long, filt.tiss, jexperiment="array", blueend = jblueend, blackend = jblackend, min.n = jmin.n, max.n = jmax.n, jmin.col = "red", jmax.col = "green", jscale=FALSE))
  if (all(fits$model == "Liver")){
    print(PlotHeatmapNconds(fits, dat.long, filt.tiss, jexperiment="array", blueend = -0.15, blackend = 0.15, min.n = -2.5, max.n = 2.5, jmin.col = "red", jmax.col = "green", jscale=FALSE))
  }
  if (min(fits$n.rhyth) > 7){
    print(PlotHeatmapNconds(subset(fits, amp.avg > 0.4), dat.long, filt.tiss, jexperiment="array", blueend = -0.15, blackend = 0.15, min.n = jmin.n, max.n = jmax.n, jmin.col = "red", jmax.col = "green", jscale=FALSE))
  }
  # plot SVD
  comp <- 1; dotsize <- 5
  genes <- as.character(unique(fits$gene))
  s <- SvdOnComplex(subset(dat.complex, gene %in% genes), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, jsize = 16)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  print(eigens$u.plot + ylab("ZT") + ggtitle(""))
  print(eigens$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))
  
  eigens.arrow2 <- GetEigens(s, period = 24, comp = comp, adj.mag = TRUE, 
                             constant.amp = 5, 
                             label.n = 5, jtitle = "", 
                             peak.to.trough = TRUE, 
                             dotsize = 1, 
                             dotshape = 18,
                             disable.text = FALSE, 
                             add.arrow = TRUE,
                             disable.repel = TRUE)
  print(eigens.arrow2$u.plot + ylab("ZT") + xlab("Log2 Fold Change") + ggtitle(""))
  print(eigens.arrow2$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))
})
dev.off()

# Heatmaps on Liver Kidney WTKO -------------------------------------------


top.models <- c("Liver_SV129,Liver_BmalKO", "Liver_SV129", "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO", "Liver_SV129,Kidney_SV129", "Kidney_SV129", "Kidney_SV129,Kidney_BmalKO")

dat.wtko$tissue <- factor(dat.wtko$tissue, levels = c("Liver", "LiverBmalKO", "Kidney", "KidneyBmalKO"))

pdf(paste0("plots/heatmaps/polarplots_heatmaps_wtko_model_selection.blacklims", jblueend, jblackend, ".minmax.", jmin.n, jmax.n, ".pdf"))
for (jmod in top.models){
  out <- PlotHeatmapNconds(subset(fits.long.filt, model == jmod), dat.wtko, filt.tiss=c(), jexperiment="RNASeq", 
                           blueend = jblueend, blackend = jblackend, min.n = jmin.n, max.n = jmax.n,
                           jmin.col = "red", jmax.col = "green", jscale=FALSE)
}
dev.off()


# Plot distribution of phases ---------------------------------------------


