# 2016-06-05
# Run Bayes Factors on Kidney and Liver from Cedric's data

rm(list=ls())

# Functions ---------------------------------------------------------------

library(dplyr)
library(ggplot2)

source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")

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


# Run bayes factors -------------------------------------------------------

gene <- "Slc44a1"
test <- MakeDesMatRunFitEnv(dat.env, gene, as.character(unique(dat$tissue)), 
                            n.rhyth.max = 2, w = 2 * pi / 24, criterion = "eb", normalize.weights = TRUE, cutoff = 1e-5, top.n = NULL, sparse = FALSE)

