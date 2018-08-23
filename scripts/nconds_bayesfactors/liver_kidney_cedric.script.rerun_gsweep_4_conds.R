# Jake Yeung
# Date of Creation: 2018-08-21
# File: ~/projects/tissue-specificity/scripts/nconds_bayesfactors/liver_kidney_cedric.script.rerun_gsweep_4_conds.R
# Rerun with 4 conditions


rm(list=ls())
setwd("/home/yeung/projects/tissue-specificity")

args <- commandArgs(trailingOnly=TRUE)
method <- args[[1]]

print(paste("METHOD:", method))

DATA="nestle"  # hogenesch|nestle
CENTER=FALSE
SCALE=FALSE
outdir <- paste0("Robjs/bayes_factors_gsweep.", Sys.Date())
dir.create(outdir)
outf <- paste0(outdir, "/fits.bayesfactors.4conds.data.", 
               DATA, ".center.", CENTER, ".scale.", SCALE, ".meth.", method, ".Robj")
print(paste("Outf:", outf))
if (file.exists(outf)) stop(paste("Outf exists, exiting", outf))


# Functions ---------------------------------------------------------------

require(dplyr)
require(ggplot2)
# require(BayesFactor)  # for checking

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")


# Load and wrangle data ---------------------------------------------------

if (DATA=="nestle"){
  load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
  dat <- dat.long; rm(dat.long)
  dat$tissue <- paste(dat$tissue, dat$geno, sep = "_")
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
n.tissues <- length(tissues.uniq)
# Fit models --------------------------------------------------------------

dat.env <- DatLongToEnvironment(dat)

# do for real
start <- Sys.time()

  print(paste("Outf:", outf))
  fits.all <- lapply(ls(dat.env), function(gene){
    MakeDesMatRunFitEnv(dat.env, gene, tissues.uniq, 
                        n.rhyth.max = n.tissues, w = 2 * pi / 24, 
                        criterion = method, normalize.weights = TRUE, 
                        cutoff = 1e-5, top.n = NULL, sparse = FALSE)
  })
  print(Sys.time() - start)
  
  fits.all.long <- lapply(fits.all, function(x){
    gene <- x$gene
    x$gene <- NULL
    fits.temp.long <- ListToLong(x, gene, top.n = 52)
  })
  fits.all.long <- do.call(rbind, fits.all.long)
  save(fits.all, fits.all.long, file=outf)

print(Sys.time() - start)

