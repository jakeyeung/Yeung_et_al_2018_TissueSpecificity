# 2015-09-28
# Use BIC but on limited set of models
# Jake Yeung

library(dplyr)
library(parallel)
setwd("~/projects/tissue-specificity")


# Source ------------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")

# Functions ---------------------------------------------------------------



GetRhythmicFormula <- function(exprs, time, intercepts, w = 2 * pi / 24, with.intercept = FALSE){
  # Get formula for full rhythmic model
  # tissue: genotype, or tissues
  # x: time
  # y: expression
  # experiment: column name of experiment (can be blank I think "")
  # intercept: if you want tissue and experiments to be your intercepts, keep it False
  
  # get rhythmic terms
  sinterm <- paste0(tissue, " : ", "sin(", w, " * ", time, ")")
  costerm <- paste0(tissue, " : ", "cos(", w, " * ", time, ")")
  rhyth.terms <- paste(sinterm, costerm, sep = " + ")
  
  if (with.intercept){
    intercept <- "1"
  } else {
    intercept <- "0"
  }
  
  intercepts <- paste(intercepts, collapse = "+")
  
  form <- as.formula(paste0(exprs, "~", paste(intercept, tissue, experiment, rhyth.terms, sep = "+")))
  return(form)
}

SubsetFullDesignMatrix <- function(des.mat, rhythmic.tissues){
  # Remake design matrix to include only "rhythmic.tissues" given the full des.mat from GetRhythmicFormula
  
  # keep intercept terms
  cols.int <- !grepl(":sin|:cos", colnames(des.mat))  # all terms without :sin or :cos are intercept terms
  
  if (length(rhythmic.tissues) == 0){
    # no rhytmic tissues
    return(des.mat[, cols.int])
  }
  
  # keep subset of rhythmic terms, depending on rhythmic tissues
  grep.str <- paste0(rhythmic.tissues, ":")  # greps sin and cos terms
  grep.str <- paste0(grep.str, collapse = "|")  # greps OR
  
  cols.rhyth <- grepl(grep.str, colnames(des.mat))
  # take intercept and subset of rhythmics as subset
  des.mat <- des.mat[, cols.int | cols.rhyth]
  return(des.mat)
}

GetDesignMatrices <- function(dat, formula){
  # Return all possible design matrices from dat
}

GetAllCombos <- function(tissues, ignore.full = TRUE){
  # get subsets of tissues possible
  # ignore.full: does not consider the full model as a combination
  tissue.combos.lst <- list()
  model.i <- 1

  if (ignore.full){
    max.rhyth <- length(tissues) - 1  # ignore full model
  } else {
    max.rhyth <- length(tissues)
  }
  
  for(i in 0:max.rhyth){
    tissues.combo <- combn(tissues, i)
    for(j in seq(ncol(tissues.combo))){
      tissue.combos.lst[[model.i]] <- tissues.combo[, j]
      model.i <- model.i + 1
    }
  }
  return(tissue.combos.lst)
}

BICWeight <- function(bic.vec){
  # Vector of BIC values, convert to weights "probability"
  # BIC = -2 log (BayesFactor)
  # therefore, BICw = exp(-0.5 * BIC) normalized across all BICs
  bic.vec.weight <- exp(-0.5 * bic.vec)
  bic.vec.weight.frac <- bic.vec.weight / sum(bic.vec.weight)
}

FitCombinations <- function(dat.gene, tiss.combos, N=3, n.cores=30){
  # Return dataframe of combinations, fit, BIC for each possible model (2^N)
  # tiss.combos: list of tissue combinations to run lienar model
  # N: return only a subset of BIC models (top 3 by default)
  tissues <- unique(dat.gene$tissue)
  gene <- dat.gene$gene[[1]]
  n.models <- 2 ^ length(tissues)
  
  fits.lst <- list()
  bics.lst <- vector(length = n.models)
  
  #   print(paste("Number of models to fit:", n.models))
  
  # BEGIN: FULL DESIGN MATRIX
  form <- GetRhythmicFormula("exprs", "time", c("tissue", "experiment"), with.intercept = FALSE)
  des.mat.full <- model.matrix(form, dat.gene)
  # END: FULL DESIGN MATRIX
  
  # BEGIN: INIT FITTING WITH FULL MODEL
  # init with fitting full model, then just update with new formula
  des.mat.sub <- des.mat.full  # rename to keep fit output names consistent
  fit.full <- lm(exprs ~ 0 + des.mat.sub, data = dat.gene)
  # END: INIT FITTING
  
  # BEGIN: FIT DIFFERENT COMBOS
  
  # PARALLEL
  #   out.lst <- mclapply(tiss.combos, function(tiss.combo){
  #     des.mat.sub <- SubsetFullDesignMatrix(des.mat.full, tiss.combo)
  #     fit.long <- dat.gene %>% do(mod = lm(exprs ~ 0 + des.mat.sub, data = .)) %>%
  #       mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","), 
  #              bic = BIC(mod[[1]]))
  #     return(fit.long)
  #   }, mc.cores = n.cores)
  
  # SERIAL
  #   out.lst <- lapply(tiss.combos, function(tiss.combo){
  #     des.mat.sub <- SubsetFullDesignMatrix(des.mat.full, tiss.combo)
  #     fit.long <- dat.gene %>% do(mod = lm(exprs ~ 0 + des.mat.sub, data = .)) %>%
  #       mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","), 
  #              bic = BIC(mod[[1]]))
  #     return(fit.long)
  #   })
  
  out.lst <- lapply(tiss.combos, function(tiss.combo){
    des.mat.sub <- SubsetFullDesignMatrix(des.mat.full, tiss.combo)
    fit.long <- dat.gene %>% 
      do(mod = FitRhythmicDesMat(., des.mat.sub)) %>%
      mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","))
    return(fit.long)
  })
  
  # END: FIT DIFFERENT COMBOS
  out.df <- do.call(rbind, out.lst)
  # out.df$bicweight <- BICWeight(out.df$bic)
  out.df$bicweight <- BICWeight(sapply(out.df$mod, function(x) x$bic))
  
  # get top 3
  out.df <- out.df[order(out.df$bicweight, decreasing = TRUE), ][1:N, ]
  out.df$gene <- gene
  gc()
  return(out.df)
}

FitRhythmicDesMat <- function(dat, des.mat){
  # Fit rhythmic with design matix
  mod <- lm(exprs ~ 0 + des.mat, data = dat)
  bic <- BIC(mod)
  return(list(fit = coef(mod), bic = bic))
}


# Load and plot on 12 conditions -------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)


# Create test set ---------------------------------------------------------

dat.sub <- subset(dat.long, gene == "Nr1d1")
clockgenes <- GetClockGenes()

tissues <- as.character(unique(dat.sub$tissue))

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

form <- GetRhythmicFormula("exprs", "time", c("tissue", "experiment"), with.intercept = FALSE)
des.mat.full <- model.matrix(form, dat.sub)

start <- Sys.time()
tiss.combos <- GetAllCombos(tiss.test, ignore.full = FALSE)

# dat.test$gene <- factor(as.character(dat.test$gene), levels = unique(dat.test$gene))
# dat.split <- split(dat.test, dat.test$gene)
dat.split <- split(dat.long, dat.long$gene)
fits.all <- mclapply(dat.split, FitCombinations, tiss.combos = tiss.combos, mc.cores = 50)
# fits.all <- lapply(dat.split, FitCombinations, tiss.combos = tiss.combos)

fits.all <- do.call(rbind, fits.all)

# fits.all <- dat.test %>%
#   group_by(gene) %>%
#   do(fits = FitCombinations(., tiss.combos, n.cores = 15))

# size of object
fits.all %>% object.size %>% print(unit = "MB")

save(fits.all, file = "Robjs/nconds.rnaseq_microarray.independent.Robj")
print(Sys.time() - start)




