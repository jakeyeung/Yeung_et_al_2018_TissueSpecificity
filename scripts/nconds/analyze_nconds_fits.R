# 2015-10-07
# Jake Yeung
# Analyze fits from vitalit

library(dplyr)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/AppendListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")

load("Robjs/dat.long.fixed_rik_genes.Robj")

# Functions ---------------------------------------------------------------

LoadFitsFromModels <- function(fitdir){
  fits.list <- expandingList()
  fit.files <- list.files(fitdir)
  for (i in seq(length(fit.files))){
    load(file.path(fitdir, fit.files[[i]]))  # fits
    fits.list$add(fits)
  }
  return(fits.list$as.list())
}

AddDf <- function(obj){
  # add arbitrary object to last column of data frame
  return(obj)
}

FitToMat <- function(fit, genename, weight, period = 24){
  # Turn fit into a long data frame
  model.name <- CoefToModelName(fit)
  # vector to long
  params <- CoefToParams(fit)
  dat.out <- data.frame(gene = genename, model = model.name, weight = weight)  # init
  dat.params <- dat.out %>%
    group_by(gene, model, weight) %>%
    do(param.list = AddDf(params))
  return(dat.params)
}

ListToLong <- function(fits.list, genename, top.n, period = 24){
  # form list of fits (sorted by best to worst) into a long dataframe
  # should be easily used to append or edit a data frame
  mats <- list()
  for (i in seq(top.n)){
    mats[[i]] <- FitToMat(fits.list[[i]]$fit, genename, fits.list[[i]]$weight.norm, period)
  }
  return(do.call(rbind, mats))
}

GetTopNModelsFromFits <- function(fits.list, top.n){
  fits.list <- NormalizeWeightsFromFits(fits.list)  # handles pesky fit.weights.sum  
  top.n.weights <- vector("numeric", length = top.n)
  top.n.tuples <- matrix(data = 0, nrow = top.n, ncol = 2)  # 2 because i and j indices
  for (i in seq(length(fits.list))){
    for (j in seq(length(fits.list[[i]]))){  # should be list of n.top or more
      weight.j <- fits.list[[i]][[j]]$weight.norm
      min.indx <- which.min(top.n.weights)
      if (length(weight.j) == 0) next
        if (top.n.weights[min.indx] < weight.j){
          # update weights and indcies
          top.n.weights[min.indx] <- weight.j
          # print(top.n.tuples[min.indx, c(1, 2)])
          top.n.tuples[min.indx, c(1, 2)] <- c(i, j)
          # print(top.n.tuples[min.indx, c(1, 2)])
        }
    }
  }
  # order top.n.tuples by top.n.weights
  top.n.tuples <- top.n.tuples[order(top.n.weights, decreasing = TRUE), ]
  
  # get the top.n tuples
  fits.top <- expandingList()
  # print(top.n.tuples)
  apply(top.n.tuples, MARGIN = 1, function(tup){
    # print(tup)
    fits.top$add(fits.list[[tup[1]]][[tup[2]]])
  })
  return(fits.top$as.list())
}

NormalizeWeightsFromFits <- function(fits.list){
  fit.weight.sum <- SumWeights(fits.list)
  for (i in seq(length(fits.list))){
    fits.list[[i]][["fit.weights.sum"]] <- NULL
    for (j in seq(length(fits.list[[i]]))){
      fits.list[[i]][[j]][["weight.norm"]] <- fits.list[[i]][[j]][["weight"]] / fit.weight.sum
    }
  }
  return(fits.list)
}


SumWeights <- function(fits.list){
  # Given list of fits, divide each weight by fit.weights.sum
  fit.weight.sum <- 0
  for (fit in fits.list){
    if (!is.null(fit$fit.weights.sum)){
      fit.weight.sum <- fit.weight.sum + fit$fit.weights.sum
    }
  }
  return(fit.weight.sum)
}

GetNrhythFromModel <- function(model){
  # from model, return number of tissues that are rhythmic 
  # e.g.: Mus;Adr,Kidney,Liver;Aorta,BFAT,Heart,Lung = 8 rhythmic tissues
  # strategy: replace ";" with ',', then get length of ',' split string
  n.rhyth <- gsub(pattern = ";", replacement = ",", x = model)
  return(length(strsplit(n.rhyth, ",")[[1]]))
}

# Load --------------------------------------------------------------------

# genedir <- "/home/yeung/projects/tissue-specificity/data/nconds2/fits_11_tiss_chunks.vitalit"
genedir <- "/home/yeung/projects/tissue-specificity/results/nconds/fits_11_tiss_3_max_cluster"

fits.long.list <- expandingList()
for (gene in list.files(genedir)){
  fitdir <- file.path(genedir, gene)
  fits <- LoadFitsFromModels(fitdir)
  fits.top <- GetTopNModelsFromFits(fits, 3)
  fits.long.gene <- ListToLong(fits.top, genename = gene, top.n = 3, period = 24)
  fits.long.list$add(fits.long.gene)
}
fits.long.list <- fits.long.list$as.list()
fits.long <- do.call(rbind, fits.long.list)
head(fits.long)

fits.long$n.rhyth <- sapply(fits.long$model, GetNrhythFromModel)

outpath <- "/home/yeung/projects/tissue-specificity/results/nconds/fits_11_tiss_max_3.Robj"
save(fits.long, file = outpath)
