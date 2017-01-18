rm(list=ls())

remove.wfat <- TRUE

library(ggplot2)
library(hash)
library(dplyr)
library(reshape2)
library(gplots)
library(penalizedLDA)
library(wordcloud)
# First source my functions I wrote
funcs.dir <- file.path('scripts', 'functions')
source(file.path(funcs.dir, 'SampleNameHandler.R'))  # for shortening sample names
source(file.path(funcs.dir, 'PcaPlotFunctions.R'))  # for visualizing PCA, periodoigrams
source(file.path(funcs.dir, 'FourierFunctions.R'))  # for periodoigrams
source(file.path(funcs.dir, 'GetTissueSpecificMatrix.R'))  # as name says
source(file.path(funcs.dir, "GrepRikGenes.R"))
source(file.path(funcs.dir, "GetTissueTimes.R"))
source(file.path(funcs.dir, "PCAFunctions.R"))
source(file.path(funcs.dir, "LoadAndHandleData.R"))
source(file.path(funcs.dir, "FitRhythmic.R"))
source(file.path(funcs.dir, "PlotGeneAcrossTissues.R"))
source(file.path(funcs.dir, "LoadArray.R"))
source(file.path(funcs.dir, "VarianceFunctions.R"))
source(file.path(funcs.dir, "FitRhythmicAcrossPeriods.R"))
source(file.path(funcs.dir, "GetClockGenes.R"))
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/OuterComplex.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/RemoveCommasBraces.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/FourierFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/LongToMat.R")


load(file = "Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/dat.fit.Robj", v=T); dat.fit.24 <- dat.fit

library(hash)

fits.best$n.rhyth.fac <- as.factor(sapply(as.numeric(fits.best$n.rhyth), function(n) NrhythToStr(n)))

filt.tiss <- c("WFAT")
load("Robjs/dat.complex.fixed_rik_genes.Robj")

if (remove.wfat){
  dat.complex <- subset(dat.complex, ! tissue %in% filt.tiss)
}

fits.rhyth <- subset(fits.best, n.params > 0)
fits.rhyth$label <- apply(fits.rhyth, 1, function(row){
  cutoff <- 1
  if (row[8] > cutoff & row[6] > 0){  # amp.avg > cutoff only for n.rhyth > 1
    return(as.character(row[1]))  # return gene
  } else {
    return("")
  }
})

# count based on amp
amp.thres <- seq(from = 0, to = max(dat.fit.24$amp), by = 0.15)

fits.best$n.rhyth.lab <- sapply(fits.best$n.rhyth, function(n){
  if (n >= 8){
    return("8-11")
  } else if (n == 1){
    return("1")
  } else if (n <= 3 & n >= 2){
    return("2-3")
  } else if (n <= 7 & n >- 4){
    return("4-7")
  } else {
    print(n)
    warning("Didnt fit any if statements")
  }
})
fits.counts.by.amp <- subset(fits.best, n.rhyth > 0) %>%
  group_by(n.rhyth.lab) %>%
  do(NGenesByAmp.long(., amp.thres, labelid = "n.rhyth.lab", varid = "amp.avg", outlabel = "n.rhyth.lab"))
ggplot(fits.counts.by.amp, aes(x = 2 * amp.thres, y = n.genes, group = n.rhyth.lab, colour = as.factor(n.rhyth.lab))) + geom_line() + 
  geom_line(size = 2) + 
  theme_bw(20) +
  labs(colour = "# Rhythmic\nTissues") + 
  theme(aspect.ratio=1, 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Avg Amplitude of Rhythmic Tissues") + ylab("# Genes") + xlim(c(0.15, 6)) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000)) + 
  geom_vline(xintercept = 2.8, linetype = "dotted") + 
  scale_colour_brewer(palette = "Spectral")