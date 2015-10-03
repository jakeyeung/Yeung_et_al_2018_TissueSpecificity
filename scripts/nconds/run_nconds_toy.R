# 2015-09-28
# Use BIC but on limited set of models
# Jake Yeung

w <- 2 * pi / 24

library(dplyr)
library(parallel)
library(hash)

setwd("~/projects/tissue-specificity")


# Source ------------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")

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

tiss.combos <- GetAllCombos(tiss.test, ignore.full = FALSE)
print(paste("Tissue combos:", length(tiss.combos)))


# Iterate n rhythmic params ------------------------------------------------------



tiss.combos <- GetAllCombos(tiss.test, ignore.full = FALSE)  # to know all combinations
my_mat.queue <- new.queue()

# init with flat model
des.mat.flat <- GetFlatModel(dat.gene)
# get rhythmic parameters which will be used for adding later: hash structure has fast lookup
des.mat.sinhash <- GetSinCombos(dat.gene, w, tiss.test, tiss.combos)
des.mat.coshash <- GetCosCombos(dat.gene, w, tiss.test, tiss.combos)

rhyth.tiss <- list(character(0))  # needs to track shared and independent parameters, e.g.: c("Liver,Kidney", "Adr") no duplicates allowed
# n.rhyth <- NRhythmicFromString(rhyth.tiss)  # number of independent rhythmic parameters perhaps? 
n.rhyth <- NRhythmicFromVector(rhyth.tiss)  # do that later? naw faster if we do it now
complement <- FilterCombos(tiss.combos, rhyth.tiss)  # large memory usage let's try not to optimize unless we need to 
# n.param <- 1  # just flat
des.mat.list <- list(mat=des.mat.flat, rhyth.tiss=rhyth.tiss, n.rhyth=n.rhyth, complement = complement)

enqueue(my_mat.queue, des.mat.list)

n.mat.submitted <- 1

# generate matrix by adding combinations of columns and adding
# those matrices into the queue
while (! is.empty(my_mat.queue)) {
  des.mat.list <- dequeue(my_mat.queue)
  # determine tissue combinations that need to be added based on rhyth.tiss
  # e.g., no need to add Liver twice, they can't have two rhythmic paramters
  
  #   tiss.combos.sub <- FilterCombos(tiss.combos, des.mat.list$rhyth.tiss)

  for (tiss.comb in des.mat.list$complement){
    # add column for each tissue combination
    tiss.key <- paste(tiss.comb, collapse = ",")
    col.new <- AddRhythmicColumns(des.mat.sinhash, des.mat.coshash, tiss.key)
    
    # append tiss.key to rhyth.tiss
    rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.comb)

    # further remove complement after having
    tiss.complement.new <- FilterCombos(des.mat.list$complement, tiss.comb)
    
    # add meta data: makes finding models easier
    n.rhyth <- des.mat.list$n.rhyth + length(tiss.comb)
    
    # make new matrix, put it into queue
    mat.new <- cbind(des.mat.list$mat, col.new)
    des.mat.list.new <- list(mat=mat.new, rhyth.tiss = rhyth.tiss, n.rhyth=n.rhyth, complement = tiss.complement.new)
    enqueue(my_mat.queue, des.mat.list.new) 
    n.mat.submitted <- n.mat.submitted + 1
  }  
} 
print(n.mat.submitted)


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




