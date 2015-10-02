# 2015-09-28
# Use BIC but on limited set of models
# Jake Yeung

library(dplyr)
library(parallel)
setwd("~/projects/tissue-specificity")


# Source ------------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")
source("scripts/functions/NcondsFunctions.R")

# Functions ---------------------------------------------------------------

GetRhythmicSharedParams2 <- function(dat.gene, rhyth.tiss, tissue.i = 2){
  # Given dat.gene and tissues that are rhythmic, return design matrices
  # for corresponding model
  w <- 2 * pi / 24
  dat.gene$is.rhyth <- apply(dat.gene, 1, function(dat.row) as.numeric(dat.row[tissue.i] %in% rhyth.tiss))
  my_mat.shared <- model.matrix(exprs ~ 0 + tissue + experiment + I(is.rhyth*(sin(w * time))), data = dat.gene)
}

ConcatenateRhythmicTissues <- function(des.mat, des.mat.rhyth.sin, des.mat.rhyth.cos, rhyth.tiss){
  # Concatenate rhythmic tissues so they have shared parameters
  grepl.str <- paste(rhyth.tiss, collapse = "|")
  des.mat.cols <- lapply(list(des.mat.rhyth.sin, des.mat.rhyth.cos), function(d.mat){
    d.mat.sub <- d.mat[, colnames(d.mat)[grepl(grepl.str, colnames(d.mat))]]
    d.mat.sub <- as.matrix(apply(d.mat.sub, 1, sum))
  })
  des.mat.cols <- do.call(cbind, des.mat.cols)
  des.mat <- cbind(des.mat, des.mat.cols)
  return(des.mat)
}

enqueue <- function(queue, new.element){
  temp.queue <- new.env()
  temp.queue$value <- new.element
  temp.queue$next.element <- queue$next.element
  queue$next.element <- temp.queue
  queue$value <- NULL
}

dequeue <- function(queue){
  value <- queue$next.element$value
  queue$next.element <- (queue$next.element)$next.element
  queue$value <- NULL
  return(value)
}

system.time({
  x <- new.env()
  for(i in 1:N){
    enqueue(x,i)
  } 
}) 

# Load and plot on 12 conditions -------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)

# remove WFAT
dat.long <- subset(dat.long, tissue != "WFAT")


# Create test set ---------------------------------------------------------

dat.sub <- subset(dat.long, gene == "Nr1d1")
clockgenes <- GetClockGenes()

tissues <- as.character(unique(dat.long$tissue))

tiss.test <- c("Liver", "Kidney", "Adr", "BFAT", "Mus")
# tiss.test <- tissues
dat.test <- subset(dat.long, gene %in% clockgenes & tissue %in% tiss.test)
# dat.test <- subset(dat.long, tissue %in% tiss.test)

# Label colnames ----------------------------------------------------------

experiment <- "experiment"
time <- "time"
exprs <- "exprs"
tissue <- "tissue"


# Get design matrices -----------------------------------------------------

# form <- GetRhythmicFormula("exprs", "time", c("tissue", "experiment"), with.intercept = FALSE)
# des.mat.full <- model.matrix(form, dat.sub)

start <- Sys.time()
tiss.combos <- GetAllCombos(tiss.test, ignore.full = FALSE)
print(paste("Tissue combos:", length(tiss.combos)))


# Iterate n rhythmic params ------------------------------------------------------

my_mat <- list()
model <- 1
models.done <- c()
des.mat.flat <- GetFlatModel(dat.gene)
my_mat$model[[model]] <- des.mat.flat
my_mat$rhyth.tiss <- ""

n.rhyth.param <- 1  # test
dat.gene <- subset(dat.long, gene == "Nr1d1")

for (tiss.comb in tiss.combos){
  # iterate each tiss combo to have n.rhyth.parameters
  
}



# dat.test$gene <- factor(as.character(dat.test$gene), levels = unique(dat.test$gene))
# dat.split <- split(dat.test, dat.test$gene)
dat.split <- split(dat.long, dat.long$gene)
fits.all <- mclapply(dat.split, FitCombinations, tiss.combos = tiss.combos, mc.cores = 50)
# start <- Sys.time()
# fits.all <- lapply(dat.split, FitCombinations, tiss.combos = tiss.combos)
fits.all <- do.call(rbind, fits.all)
# print(Sys.time() - start)

# fits.all <- dat.test %>%
#   group_by(gene) %>%
#   do(fits = FitCombinations(., tiss.combos, n.cores = 15))

# size of object
fits.all %>% object.size %>% print(unit = "MB")

save(fits.all, file = "Robjs/nconds.rnaseq_microarray.independent.Robj")
print(Sys.time() - start)




