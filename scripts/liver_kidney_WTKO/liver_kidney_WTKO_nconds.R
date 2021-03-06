# 2016-06-23
# Run nconds on 4 conditions: Liver Kidney WT and KO

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")


require(dplyr)
require(ggplot2)
require(BayesFactor)  # for checking
require(numbers)  # for Bell number

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LiverKidneyFunctions.R")

args <- commandArgs(trailingOnly=TRUE)
method <- args[[1]]
print(paste("METHOD:", method))
CENTER=FALSE
SCALE=FALSE
removesamps=TRUE
stagger=TRUE
# outf <- paste0("Robjs/liver_kidney_atger_nestle/nconds/fits.nconds.center.", CENTER, ".scale.", SCALE, ".meth.", method, ".Robj")
outf <- paste0("Robjs/liver_kidney_atger_nestle/nconds/fits.nconds.stagger.", stagger, ".removesamps.", removesamps, ".meth.", method, ".Robj")

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)


# Filter bad genes --------------------------------------------------------

# Can filter bad genes now, but better to do it downstream!
# dat.long <- RemoveLowExprsPseudoShortGenes(dat.long, gbiotype = "protein_coding", glength = 250, show.plot = FALSE)

# Prepare data to fit models ----------------------------------------------

dat.long <- CollapseTissueGeno(dat.long)
if (removesamps){
  if (stagger){
    # Stagger removes ZT2 from Kidney, which is an outlier when we use Kallisto
    dat.long <- StaggeredTimepointsLivKid(dat.long)
  } else {
    dat.long <- SameTimepointsLivKid(dat.long)
  }
}

dat.long <- subset(dat.long, select = c(gene, tissue, time, experiment, exprs))

tissues.uniq <- unique(as.character(dat.long$tissue))
dat.env <- DatLongToEnvironment(dat.long)


# Run ---------------------------------------------------------------------

start <- Sys.time()

print(paste("Outf:", outf))
fits.all <- lapply(ls(dat.env), function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, tissues.uniq, 
                      n.rhyth.max = length(tissues.uniq), w = 2 * pi / 24, 
                      criterion = method, normalize.weights = TRUE, 
                      cutoff = 1e-5, top.n = NULL, sparse = FALSE)
})
print(Sys.time() - start)

n.combos <- bell(length(tissues.uniq) + 1)
fits.all.long <- lapply(fits.all, function(x){
  gene <- x$gene
  x$gene <- NULL
  fits.temp.long <- ListToLong(x, gene, top.n = n.combos)
})
fits.all.long <- do.call(rbind, fits.all.long)
save(fits.all.long, file=outf)

print(Sys.time() - start)

