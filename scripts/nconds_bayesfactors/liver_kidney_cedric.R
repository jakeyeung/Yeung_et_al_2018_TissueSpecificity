# 2016-06-05
# Run Bayes Factors on Kidney and Liver from Cedric's data

rm(list=ls())

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

TissueFromCname <- function(cname){
  # "X0_1 -> Kidney"
  # "X0_1_liver -> Liver"
  jsplit <- strsplit(cname, "_")[[1]]
  if (jsplit[[length(jsplit)]] == "liver"){
    tiss <- "Liver"
  } else {
    tiss <- "Kidney"
  }
  return(tiss)
} 

TimeFromCname <- function(cname, rm.first.char=TRUE){
  # X0_1 -> 0 + (1 - 1) * 24
  # X22_2 -> 22 + (2 - 1) * 24
  if (rm.first.char){
    cname <- substr(cname, start = 2, stop = nchar(cname))
  }
  time.base <- as.numeric(strsplit(cname, "_")[[1]][[1]])
  time.multiplier <- as.numeric(strsplit(cname, "_")[[1]][[2]]) - 1
  time <- time.base + time.multiplier * 24
  return(time)
}


# Load and wrangle data ---------------------------------------------------


inf <- "/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/DAT_l_l.txt"
dat.mat <- read.table(inf, header = TRUE)
genes <- dat.mat$name
exprs <- subset(dat.mat, select = c(-name))
tissues <- rep(sapply(colnames(exprs),TissueFromCname, USE.NAMES = FALSE), each = nrow(exprs))
times <- rep(sapply(colnames(exprs), TimeFromCname, USE.NAMES = FALSE), each = nrow(exprs))

# make long
dat <- data.frame(gene = dat.mat$name, 
                  tissue = tissues,
                  time = times,
                  exprs = unlist(exprs),
                  experiment = "rnaseq")

# Fit models --------------------------------------------------------------

dat.env <- DatLongToEnvironment(dat)


# Plot genes --------------------------------------------------------------

PlotGeneAcrossTissues(subset(dat, gene == "Lrriq3"))


# Run bayes factors -------------------------------------------------------
tissues <- as.character(unique(dat$tissue))

# test
jgene <- "Slc44a1"
jgene <- "Lrriq3"
jgene <- "Cdc20"
jgene <- ls(dat.env)[[10000]]
jgene <- "Ppp3ca"
jgene <- "Alb"
jgene <- "Slc44a1"
test <- MakeDesMatRunFitEnv(dat.env, jgene, tissues, 
                            n.rhyth.max = 2, w = 2 * pi / 24, criterion = "zf", normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = FALSE)
PlotGeneAcrossTissues(subset(dat, gene == jgene))

exp(linearReg.R2stat(N=48, p=2, R2=0.999199925619844, rscale = 0.353553390593274)[['bf']])
exp(linearReg.R2stat(N=48, p=2, R2=0.8, rscale = 1)[['bf']])
GetBayesFactor(N = 48, p = 2, rsquared = 0.8, method = "zf")

# do for real
start <- Sys.time()
# Rprof()

# genes <- c("Lrriq3", "Cdc20")
fits.all <- lapply(ls(dat.env), function(gene){
  print(gene)
  # fits.all <- mclapply(genes, function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, tissues, 
                      n.rhyth.max = 2, w = 2 * pi / 24, 
                      criterion = "zf", normalize.weights = TRUE, 
                      cutoff = 1e-5, top.n = NULL, sparse = FALSE)
})
print(Sys.time() - start)


# genes <- c("Lrriq3", "Cdc20")
fits.all <- mclapply(ls(dat.env), function(gene){
# fits.all <- mclapply(genes, function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, tissues, 
                      n.rhyth.max = 2, w = 2 * pi / 24, 
                      criterion = "zf", normalize.weights = TRUE, 
                      cutoff = 1e-5, top.n = NULL, sparse = FALSE)
}, mc.cores = 5)
print(Sys.time() - start)

fits.all.long <- lapply(fits.all, function(x){
  gene <- x$gene
  x$gene <- NULL
  fits.temp.long <- ListToLong(x, gene, top.n = 5)
})
fits.all.long <- do.call(rbind, fits.all.long)
print(Sys.time() - start)

