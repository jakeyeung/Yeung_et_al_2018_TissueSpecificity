# 2016-06-16
# Choose optimal g by finding g with largest number of Kidney,Liver genes.

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FitRhythmicAcrossPeriods.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/ListFunctions.R")

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
# Load --------------------------------------------------------------------

dat <- LoadLivKid()
foutpath <- "Robjs/liver_kidney/dat.fit.periods.genome_wide.Robj"
load(foutpath, verbose = T)

indir <- "/home/yeung/projects/tissue-specificity/Robjs/bayes_factors_gsweep"
inf <- list.files(indir)
fits.long <- expandingList()
for (f in inf){
  print(f)
  method <- strsplit(f, split = "\\.")[[1]][[11]]
  fpath <- file.path(indir, f)
  load(fpath)
  fits.all.long$method <- method
  fits.long$add(fits.all.long)
}
fits.long <- fits.long$as.list()
fits.long <- do.call(rbind, fits.long)
fits.long$g <- as.numeric(sapply(fits.long$method, function(m) strsplit(m, "=")[[1]][[2]]))

# Filter ------------------------------------------------------------------

cutoff <- -0.187448

dat.mean <- dat %>%
  group_by(gene) %>%
  summarise(exprs.mean = mean(exprs))

dat.mean.bytiss <- dat %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs))

genes.filt <- as.character(subset(dat.mean, exprs.mean < cutoff)$gene)

dat <- subset(dat, ! gene %in% genes.filt)

# get best model
fits.long.filt <- subset(fits.long, ! gene %in% genes.filt) %>%
  group_by(gene, g) %>%
  filter(weight.raw == min(weight.raw))

# fourier
omega <- 2 * pi / 24
dat.freq <- dat %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))


# Count models ------------------------------------------------------------

fits.count <- fits.long.filt %>%
  group_by(g, model) %>%
  summarize(n.genes = length(model))



# Count by temporal variance ----------------------------------------------

library(hash)

dat.freq.tvar <- dat.freq %>%
  group_by(gene) %>%
  summarize(tvar = sum(Mod(exprs.transformed * 2) ^ 2))
temp.var <- hash(as.character(dat.freq.tvar$gene), dat.freq.tvar$tvar)

fits.long.filt$tvar <- sapply(as.character(fits.long.filt$gene), function(g) temp.var[[g]])

fits.count <- fits.long.filt %>%
  group_by(g, model) %>%
  summarize(n.genes = length(model),
            tvar = sum(tvar))

tvar.flat <- hash(as.character(subset(fits.count, model == "")$g), subset(fits.count, model == "")$tvar)
fits.count$tvar.flat <- sapply(as.character(fits.count$g), function(g) tvar.flat[[g]])

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")
pdf("plots/liver_kidney_math710/gsweep.pdf")
ggplot(fits.count, aes(x = g, y = n.genes, colour = model, group = model)) + geom_point(size = 3) + geom_line() + xlim(0, 5001) +
  theme_bw() + 
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  ylab("# genes")
  

ggplot(fits.count, aes(x = g, y = tvar, colour = model, group = model)) + geom_point(size = 3) + geom_line() + xlim(0, 5001) +
  theme_bw() + 
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  scale_colour_manual(values=cbPalette) +
  ylab("24h Spectral Power") + 
  geom_vline(xintercept = 1000)

ggplot(fits.count, aes(x = g, y = tvar, colour = model, group = model)) + geom_point(size = 3) + geom_line() + xlim(0, 10001) +
  theme_bw() + 
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  scale_colour_manual(values=cbPalette) +
  ylab("24h Spectral Power") + 
  geom_vline(xintercept = 1000)
  
# ratio of flat model
ggplot(fits.count, aes(x = g, y = tvar / tvar.flat, colour = model, group = model)) + geom_point(size = 3) + geom_line() + xlim(0, 1001) +
  theme_bw() + 
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  scale_colour_manual(values=cbPalette) +
  ylab("tvar / tvarflat")
dev.off()