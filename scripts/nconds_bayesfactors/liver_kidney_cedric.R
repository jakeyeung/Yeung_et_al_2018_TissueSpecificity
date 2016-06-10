# 2016-06-05
# Run Bayes Factors on Kidney and Liver from Cedric's data

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

# Functions ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(BayesFactor)  # for checking

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")


# Load and wrangle data ---------------------------------------------------

dat <- LoadLivKid()

# Fit models --------------------------------------------------------------

dat.env <- DatLongToEnvironment(dat)


# Plot genes --------------------------------------------------------------

# PlotGeneAcrossTissues(subset(dat, gene == "Lrriq3"))


# Run bayes factors -------------------------------------------------------
# tissues <- as.character(unique(dat$tissue))
# 
# # test
# jgene <- "Slc44a1"
# jgene <- "Lrriq3"
# jgene <- "Cdc20"
# jgene <- ls(dat.env)[[10000]]
# jgene <- "Ppp3ca"
# jgene <- "Alb"
# jgene <- "Slc44a1"
# jgene <- "4930547N16Rik"
# test <- MakeDesMatRunFitEnv(dat.env, jgene, tissues.uniq, n.rhyth.max = 2, w = 2 * pi / 24, 
#                             criterion = "hyperg", normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = FALSE)
# PlotGeneAcrossTissues(subset(dat, gene == jgene))

# do for real
start <- Sys.time()

for (method in c("zf", "eb", "hyperg", "BIC", "AIC")){
# for (method in c("AIC")){
  print(paste("method:", method))
  outf <- paste0("Robjs/fits.bayesfactors.livkid.meth.", method, ".Robj")
  print(paste("Outf:", outf))
  fits.all <- lapply(ls(dat.env), function(gene){
    MakeDesMatRunFitEnv(dat.env, gene, tissues.uniq, 
                        n.rhyth.max = 2, w = 2 * pi / 24, 
                        criterion = method, normalize.weights = TRUE, 
                        cutoff = 1e-5, top.n = NULL, sparse = FALSE)
  })
  print(Sys.time() - start)
  
  fits.all.long <- lapply(fits.all, function(x){
    gene <- x$gene
    x$gene <- NULL
    fits.temp.long <- ListToLong(x, gene, top.n = 5)
  })
  fits.all.long <- do.call(rbind, fits.all.long)
  save(fits.all.long, file=outf)
}

print(Sys.time() - start)

