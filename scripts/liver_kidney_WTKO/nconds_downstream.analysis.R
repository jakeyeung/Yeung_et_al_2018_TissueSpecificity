# 2016-06-23
# Jake Yeung

rm(list=ls())

library(dplyr)
library(ggplot2)
library(hash)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/BiomartFunctions.R")

load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)
dat.orig <- dat.long

dat.long <- CollapseTissueGeno(dat.long)
dat.long <- SameTimepointsLivKid(dat.long)

# filter NA changes
dat.long <- subset(dat.long, !is.na(gene))

# Filter lowly expressed genes --------------------------------------------

# Remove pseudogenes
genes <- unique(as.character(dat.long$gene))
genes.biotype <- AnnotatePseudogenes(genes, return.original = FALSE)
genes.length <- AnnotateTranscriptLength(genes, return.original = FALSE)
biotype.hash <- hash(genes, genes.biotype)
length.hash <- hash(genes, genes.length)
dat.long$gbiotype <- sapply(as.character(dat.long$gene), function(g) biotype.hash[[g]])
dat.long$glength <- sapply(as.character(dat.long$gene), function(g) length.hash[[g]])
dat.long <- subset(dat.long, gbiotype == "protein_coding" & glength > 250)


dat.mean <- dat.long %>%
  group_by(gene) %>%
  summarise(exprs.max = quantile(exprs, probs = 0.9))

plot(density(dat.mean$exprs.max))
jcutoff <- 1
abline(v=jcutoff)

genes.cut <- as.character(subset(dat.mean, exprs.max <= jcutoff)$gene)

dat.long <- subset(dat.long, ! gene %in% genes.cut)

# Project to Frequency ----------------------------------------------------


omega <- 2 * pi / 24
dat.freq <- dat.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))

s <- SvdOnComplex(dat.freq, value.var = "exprs.transformed")

for (i in seq(1)){
  eigens <- GetEigens(s, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
}

# All periods -------------------------------------------------------------


loadfile <- "Robjs/liver_kidney_atger_nestle/dat.complex.all_T.Robj"
if (file.exists(loadfile)){
  load(loadfile)
} else {
  periods <- rep(48, 6) / seq(1, 6)  # 48/1, 48/2 ... 48/12
  
  library(parallel)
  dat.complexes <- mclapply(periods, function(period, dat.long){
    omega <- 2 * pi / period
    dat.tmp <- dat.long %>%
      group_by(gene, tissue) %>%
      do(ProjectToFrequency2(., omega, add.tissue=TRUE))
    dat.tmp$period <- period
    return(dat.tmp)
  }, dat.long = dat.long, mc.cores = length(periods))
  
  
  dat.complex.all_T <- do.call(rbind, dat.complexes)
  save(dat.complex.all_T, file = "Robjs/liver_kidney_atger_nestle/dat.complex.all_T.Robj")
  rm(dat.complexes)
}


# By clusters -------------------------------------------------------------

jmeth <- "g=1001"
i <- 1

genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% c("Kidney_SV129,Liver_SV129"))$gene)

s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)


genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% c("Kidney_SV129"))$gene)

s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)


genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% c("Liver_SV129"))$gene)

s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = i, label.n = 30, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)



# Count by temporal variance ----------------------------------------------

fits.g <- subset(fits.long.filt, !gene %in% genes.cut & !method %in% c("zf"))
fits.g$g <- sapply(fits.g$method, function(m){
  g = tryCatch({
    g <- as.numeric(strsplit(m, "=")[[1]][[2]])
  }, error = function(e) {
    g <- m
  })
  return(g)
})

by.noise <- TRUE

if (!by.noise){
  # for 24hr variance
  dat.freq.tvar <- subset(dat.freq) %>%
    group_by(gene) %>%
    summarize(tvar = sum(Mod(exprs.transformed * 2) ^ 2))
  temp.var <- hash(as.character(dat.freq.tvar$gene), dat.freq.tvar$tvar)
  jylab <- "24h Spectral Power"
} else {
  # for 16 and 9.6 hour variance
  noise.components <- periods[which(24 %% periods != 0)]
  dat.freq.tvar <- subset(dat.complex.all_T, period %in% noise.components) %>%
    group_by(gene) %>%
    summarize(tvar = sum(Mod(exprs.transformed * 2) ^ 2))
  temp.var <- hash(as.character(dat.freq.tvar$gene), dat.freq.tvar$tvar)
  jylab <- paste0(paste(noise.components, collapse=","), "Spectral Power")
}

fits.g$tvar <- sapply(as.character(fits.g$gene), function(g){
  tvar <- temp.var[[g]]
  if (is.null(tvar)){
    return(0)
  } else {
    return(tvar)
  }
})

fits.count <- fits.g %>%
  group_by(g, model) %>%
  summarize(n.genes = length(model),
            tvar = sum(tvar))

tvar.flat <- hash(as.character(subset(fits.count, model == "")$g), subset(fits.count, model == "")$tvar)
fits.count$tvar.flat <- sapply(as.character(fits.count$g), function(g) tvar.flat[[g]])

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")

jmodels <- c("Liver_SV129", "Kidney_SV129", "Kidney_SV129,Liver_SV129", "Kidney_SV129;Liver_SV129", "")
jmodels <- c("Liver_SV129", "Kidney_SV129", "Kidney_SV129,Liver_SV129", "Kidney_SV129;Liver_SV129")

bic.var <- subset(fits.count, g == "BIC" & model %in% jmodels)

ggplot(subset(fits.count, model %in% jmodels & g != "BIC"), aes(x = as.numeric(g), y = tvar, colour = model, group = model)) + geom_point(size = 3) + geom_line() + xlim(0, 5001) + 
  theme_bw() + 
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  scale_colour_manual(values=cbPalette) +
  ylab(jylab) + 
  xlab("g (larger g favors simpler models)") + 
  geom_hline(yintercept=bic.var$tvar[[1]], colour=cbPalette[[1]], linetype = "dotted") + 
  geom_hline(yintercept=bic.var$tvar[[2]], colour=cbPalette[[2]], linetype = "dotted") + 
  geom_hline(yintercept=bic.var$tvar[[3]], colour=cbPalette[[3]], linetype = "dotted") + 
  geom_hline(yintercept=bic.var$tvar[[4]], colour=cbPalette[[4]], linetype = "dotted") + 
  geom_vline(xintercept = 1000)


