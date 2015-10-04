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
source("scripts/functions/ListFunctions.R")

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

tissues <- as.character(unique(dat.long$tissue))
# tiss.test <- c("Liver", "Kidney", "Adr", "BFAT", "Mus", "Aorta", "WFAT")
# tiss.test <- tissues
tiss.test <- c("Liver", "Kidney")
dat.gene <- subset(dat.long, gene == "Nr1d1" & tissue %in% tiss.test)
dat.gene$tissue <- factor(as.character(dat.gene$tissue), levels = tiss.test)

my_mat <- MakeRhythmicDesignMatrices(dat.gene, simplify = TRUE)

my_fits <- FitModels(dat.gene, my_mat)

fit <- lm(exprs ~ 0 + my_mat[[5]]$mat, dat.gene)
start <- Sys.time()
for (i in 1:20000){
  fit.faster <- lm.fit(x = my_mat[[1]]$mat, y = dat.gene$exprs)
}
print(Sys.time() - start)
fit.faster

tiss.combos <- GetAllCombos(tiss.test, ignore.full = FALSE)
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
des.mats <- expandingList() 
des.mats$add(des.mat.list)

# need to track models that we have done, so we eliminate "permutations" like c("Liver", "Kidney") and c("Kidney", "Liver) models
# use hash for speed
models.done <- hash()

# generate matrix by adding combinations of columns and adding
# those matrices into the queue
while (! is.empty(my_mat.queue)) {
  des.mat.list <- dequeue(my_mat.queue)
  # determine tissue combinations that need to be added based on rhyth.tiss
  # e.g., no need to add Liver twice, they can't have two rhythmic paramters
  
  for (tiss.comb in des.mat.list$complement){
    # add column for each tissue combination
    tiss.key <- paste(tiss.comb, collapse = ",")
    
    # append tiss.key to rhyth.tiss
    rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.key)  # form list("Adr,Kidney", "Mus")
    
    # check if this tissue combination has been already submitted into queue (but in different permutation)
    # track models we have done globally
    modelname <- MakeModelName(rhyth.tiss)
    if (! is.null(models.done[[modelname]])){
      # this is a permutation of an already done combo, skip
#       print(rhyth.tiss)
#       print(paste('Skipping', modelname))
      next
    }
    
    col.new <- AddRhythmicColumns(des.mat.sinhash, des.mat.coshash, tiss.key)
    
    #     rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.key)
    

    # further remove complement after having
    tiss.complement.new <- FilterCombos(des.mat.list$complement, tiss.comb)
    
    # add meta data: makes finding models easier
    n.rhyth <- des.mat.list$n.rhyth + length(tiss.comb)
    
    # make new matrix, put it into queue
    mat.new <- cbind(des.mat.list$mat, col.new)
    des.mat.list.new <- list(mat=mat.new, rhyth.tiss = rhyth.tiss, n.rhyth=n.rhyth, complement = tiss.complement.new)
    enqueue(my_mat.queue, des.mat.list.new) 
    models.done[[modelname]] <- TRUE  # we dont want to redo permutations of same models
    n.mat.submitted <- n.mat.submitted + 1
    des.mats$add(des.mat.list.new)
  }  
} 
print(n.mat.submitted)
des.mats.list <- des.mats$as.list()



