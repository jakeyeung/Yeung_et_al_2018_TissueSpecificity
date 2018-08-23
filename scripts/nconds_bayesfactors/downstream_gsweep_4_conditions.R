# Jake Yeung
# Date of Creation: 2018-08-22
# File: ~/projects/tissue-specificity/scripts/nconds_bayesfactors/downstream_gsweep_4_conditions.R
# Downstream gsweep of 4 conditions

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(hash)

source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/FitRhythmicAcrossPeriods.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/ListFunctions.R")

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")


# Load files --------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/Robjs/bayes_factors_gsweep.2018-08-22"

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

savedir <- "/home/yeung/projects/tissue-specificity/Robjs/bayes_factors_gsweep.2018-08-22/concatenated_Robj"
save(fits.long, file = file.path(savedir, "fits.long.gsweep_4_conditions.Robj"))


# Load data ---------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
dat <- dat.long; rm(dat.long)
dat$tissue <- paste(dat$tissue, dat$geno, sep = "_")

load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid

# Downstream --------------------------------------------------------------

cutoff <- -0.187448

dat.mean <- dat %>%
  group_by(gene) %>%
  summarise(exprs.mean = mean(exprs))

dat.mean.bytiss <- dat %>%
  group_by(gene, tissue) %>%
  summarise(exprs.mean = mean(exprs))

genes.filt <- c()
# genes.filt <- as.character(subset(dat.mean, exprs.mean < cutoff)$gene)
# dat <- subset(dat, ! gene %in% genes.filt)

# get best model
fits.long.filt <- subset(fits.long, ! gene %in% genes.filt) %>%
  group_by(gene, g) %>%
  filter(weight.raw == min(weight.raw))

save(fits.long, file = file.path(savedir, "fits.long.filt.gsweep_4_conditions.Robj"))


# Count models ------------------------------------------------------------

fits.count <- fits.long.filt %>%
  group_by(g, model) %>%
  summarize(n.genes = length(model))


# Temporal variance -------------------------------------------------------

dat.freq.tvar <- dat.freq %>%
  group_by(gene) %>%
  summarize(tvar = sum(Mod(exprs.transformed * 2) ^ 2))
# temp.var <- hash(as.character(dat.freq.tvar$gene), dat.freq.tvar$tvar)


fits.long.filt <- left_join(fits.long.filt, dat.freq.tvar)  # some NAs which we should remove

tvar.flat <- hash(as.character(subset(fits.count, model == "")$g), subset(fits.count, model == "")$tvar)
fits.count$tvar.flat <- sapply(as.character(fits.count$g), function(g){
  jflat <- tvar.flat[[g]]
  if (is.null(jflat)){
    jflat <- NA
  }
  return(jflat)
})


fits.count <- fits.long.filt %>%
  filter(!is.na(tvar)) %>%
  group_by(g, model) %>%
  summarize(n.genes = length(model),
            tvar = sum(tvar))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")

# jmodels <- c("Kidney_SV129,Liver_SV129", "Kidney_BmalKO,Kidney_SV129,Liver_BmalKO,Liver_SV129", "Liver_SV129", "Kidney_SV129")
jmodels <- c("Kidney_BmalKO,Kidney_SV129,Liver_BmalKO,Liver_SV129")
jmodels.wflat <- c(jmodels, "")




# Plots -------------------------------------------------------------------

pdf("/home/yeung/projects/tissue-specificity/plots/gsweep_4_conditions_test/gsweep_4_conditions_optimal_g.pdf", useDingbats = FALSE)

m.tvar <- ggplot(subset(fits.count, model %in% jmodels) , aes(x = g, y = tvar, colour = model, group = model)) + geom_point(size = 3) + geom_line() +
  theme_bw(24) +
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values=cbPalette) +
  ylab("24h Spectral Power") +
  geom_vline(xintercept = 1000)

print(m.tvar)
print(m.tvar + theme_bw() + theme(legend.position = "right"))

m.genecount <- ggplot(subset(fits.count, model %in% jmodels), aes(x = g, y = n.genes, colour = model, group = model)) + geom_point(size = 3) + geom_line() +
  theme_bw(24) +
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values=cbPalette) +
  ylab("# genes") + 
  geom_vline(xintercept = 1000)

print(m.genecount)
print(m.genecount + theme_bw() + theme(legend.position = "right"))

m.tvarnorm <- ggplot(subset(fits.count, model %in% jmodels), aes(x = g, y = tvar / tvar.flat, colour = model, group = model)) + geom_point(size = 3) + geom_line() +
  theme_bw(24) +
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom") +
  scale_colour_manual(values=cbPalette) +
  ylab("tvar / tvarflat") + 
  geom_vline(xintercept = 1000)

print(m.tvarnorm)
print(m.tvarnorm + theme_bw() + theme(legend.position = "right"))

dev.off()

