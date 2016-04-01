library(dplyr)
library(parallel)
library(ggplot2)

rm(list=ls())

load("Robjs/dat.long.fixed_rik_genes.Robj")

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")


# Subset liver to kidney --------------------------------------------------

jtiss <- c("Liver", "Kidney")

dat.merged <- subset(dat.long, tissue %in% jtiss)
dat.merged$tissue <- factor(as.character(dat.merged$tissue), levels = c("Kidney", "Liver"))
dat.merged$gene <- factor(as.character(dat.merged$gene), levels = unique(dat.merged$gene))

print("Chunking data to environment")
dat.env <- DatLongToEnvironment(dat.merged)

start <- Sys.time()

fits.all <- mclapply(ls(dat.env), function(gene){
  dat.gene <- get(gene, envir = dat.env)
  MakeDesMatRunFit(dat.gene, gene, jtiss, n.rhyth.max = 2, w = 2 * pi / 24, criterion = "BIC", normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = FALSE)
}, mc.cores = 5)
print(Sys.time() - start)

fits.all.long <- lapply(fits.all, function(x){
  gene <- x$gene
  x$gene <- NULL
  fits.temp.long <- ListToLong(x, gene, top.n = 5)
})
fits.all.long <- do.call(rbind, fits.all.long)
print(Sys.time() - start)



fits.all.long$n.params <- sapply(fits.all.long$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.all.long$n.rhyth <- sapply(fits.all.long$model, GetNrhythFromModel)
fits.all.long$amp.avg <- mapply(GetAvgAmpFromParams, fits.all.long$param.list, fits.all.long$model)
fits.all.long$phase.sd <- mapply(GetSdPhaseFromParams, fits.all.long$param.list, fits.all.long$model)
fits.all.long$phase.maxdiff <- mapply(GetMaxPhaseDiffFromParams, fits.all.long$param.list, fits.all.long$model)

save(fits.all.long, file = "Robjs/fits.liver_kidney.Robj")

