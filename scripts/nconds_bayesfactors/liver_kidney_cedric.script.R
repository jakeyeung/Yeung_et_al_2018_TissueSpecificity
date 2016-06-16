# 2016-06-05
# Run Bayes Factors on Kidney and Liver from Cedric's data

rm(list=ls())

args <- commandArgs(trailingOnly=TRUE)
method <- args[[1]]

print(paste("METHOD:", method))

DATA="nestle"  # hogenesch|nestle
CENTER=FALSE
SCALE=FALSE

setwd("/home/yeung/projects/tissue-specificity")

# Functions ---------------------------------------------------------------

require(dplyr)
require(ggplot2)
require(BayesFactor)  # for checking

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")


# Load and wrangle data ---------------------------------------------------

if (DATA=="nestle"){
  dat <- LoadLivKid()
} else if (DATA=="hogenesch"){
  load("Robjs/dat.long.fixed_rik_genes.Robj"); dat <- dat.long; rm(dat.long)
  # run on array and rnaseq
  dat <- subset(dat, tissue %in% c("Liver", "Kidney"))
  dat$tissue <- factor(as.character(dat$tissue), levels = c("Kidney", "Liver"))
  dat.mean <- subset(dat, experiment=="rnaseq") %>% 
    group_by(tissue, gene) %>%
    summarise(exprs.mean = mean(exprs))
  cutoff <- 4
  genes.exprs <- as.character(subset(dat.mean, exprs.mean >= cutoff)$gene)
  dat <- subset(dat, gene %in% genes.exprs)
  print(paste("Expressed genes:", length(unique(dat$gene))))
} else {
  print(paste("DATA:", DATA))
  stop("DATA must be nestle or hogenesch")
}

if (CENTER | SCALE){
  dat <- dat %>% 
    group_by(gene) %>%
    mutate(exprs = scale(exprs, center=CENTER, scale=SCALE))
}

print(head(dat))

tissues.uniq <- unique(as.character(dat$tissue))

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

# for (method in c("zf", "eb", "hyperg", "BIC", "AIC")){
# for (method in c("eb")){
# for (method in c("g=250", "g=500", "g=750", "g=1000")){
  print(paste("method:", method))
  outf <- paste0("Robjs/fits.bayesfactors.livkid.data.", DATA, ".center.", CENTER, ".scale.", SCALE, ".meth.", method, ".Robj")
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
# }

print(Sys.time() - start)

