# 2016-05-31
# bayes_factors.R
# Do model selection using bayes factors instead of BIC

rm(list=ls())

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)


# Do Liver vs Kidney ------------------------------------------------------

dat.sub <- subset(dat.long, tissue %in% c("Liver", "Kidney"))
dat.sub$tissue <- factor(as.character(dat.sub$tissue), levels = c("Liver", "Kidney"))

dat.env <- DatLongToEnvironment(dat.sub)
# 
# start <- Sys.time()
# # Rprof()
# fits.all <- mclapply(ls(dat.env), function(gene){
#   MakeDesMatRunFitEnv(dat.env, gene, as.character(unique(dat.sub$tissue)), n.rhyth.max = 2, w = 2 * pi / 24, criterion = "BIC", normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = FALSE)
# }, mc.cores = 5)
# print(Sys.time() - start)

# debug
gene <- "Slc44a1"
test <- MakeDesMatRunFitEnv(dat.env, gene, as.character(unique(dat.sub$tissue)), n.rhyth.max = 2, w = 2 * pi / 24, criterion = "eb", normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = FALSE)

# check
N <- 64
p <- 8
# R2 <- 0.89170133186249
R2 <- 0.10133186249
R2 <- 0.70133186249

library(BayesFactor)
# exp(linearReg.R2stat(N=64, p=8, R2=0.89170133186249, rscale = 0.353553390593274)[['bf']])
b1 <- linearReg.R2stat(N=64, p=8, R2=R2, rscale = 1)[['bf']]
b2 <- GetBayesFactor(N, p, R2, method = "zf", plot.integrand = FALSE)
print(paste(b1, b2))
print(b2 - b1)

GetBayesFactor(N, p, R2, method = "hyperg", plot.integrand = FALSE)
GetBayesFactor(N, p, R2, method = "eb", plot.integrand = FALSE)
