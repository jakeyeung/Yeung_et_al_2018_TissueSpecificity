# 2015-09-28
# Use BIC but on limited set of models
# Jake Yeung

library(dplyr)

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

FitCombinations <- function(dat.gene, tiss.combos, N=3){
  # Return dataframe of combinations, fit, BIC for each possible model (2^N)
  # tiss.combos: list of tissue combinations to run lienar model
  # N: return only a subset of BIC models (top 3 by default)
  tissues <- unique(dat.gene$tissue)
  n.models <- 2 ^ length(tissues)
  
  fits.lst <- list()
  bics.lst <- vector(length = n.models)
  
  #   print(paste("Number of models to fit:", n.models))
  
  # BEGIN: FULL DESIGN MATRIX
  form <- GetRhythmicFormula("exprs", "time", c("tissue", "experiment"), with.intercept = FALSE)
  des.mat.full <- model.matrix(form, dat.gene)
  # END: FULL DESIGN MATRIX
  
  # BEGIN: INIT FITTING WITH FULL MODEL
  model.i <- 1  # keep track of model, our row index
  # init with fitting full model, then just update with new formula
  des.mat.sub <- des.mat.full  # rename to keep fit output names consistent
  fit.full <- lm(exprs ~ 0 + des.mat.sub, data = dat.gene)
  #   fit.full.bic <- BIC(fit.full)
  #   fits.lst[[model.i]] <- fit.full
  #   bics.lst[model.i] <- fit.full.bic
  #   model.i <- model.i + 1
  # END: INIT FITTING
  
  # BEGIN: FIT DIFFERENT COMBOS
  out.lst <- lapply(tiss.combos, function(tiss.combo){
    des.mat.sub <- SubsetFullDesignMatrix(des.mat.full, tiss.combo)
    fit.long <- dat.gene %>% do(mod = lm(exprs ~ 0 + des.mat.sub, data = .)) %>%
      mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","), 
             bic = BIC(mod[[1]]))
  #                                 mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","),
  #                                        bic = BIC(mod))
  #     fit.long <- do(mod = lm(exprs ~ 0 + des.mat.sub, data = dat.gene), .data = dat.gene) %>%
  #       mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","),
  #              bic = BIC(mod))
  #     fit.long <- do(mod = update(fit.full, formula. = exprs ~ 0 + des.mat.sub), .data = des.mat.sub)
  #                      mutate(rhyth.tiss=paste0(tiss.combo, collapse = ","),
  #                             bic=BIC(mod)))
    return(fit.long)
      
  #     fit.i <- update(fit.full, formula. = exprs ~ 0 + des.mat.sub)
  #     print(summary(fit.i))
  #     fit.i.bic <- BIC(fit.i)
  #     fits.lst[[model.i]] <- fit.i
  #     bics.lst[model.i] <- fit.i.bic
  #     model.i <- model.i + 1
  #     return(list(fits=fits.lst, bics=bics.lst))
  })
  #   for (tiss.combo in tiss.combos){
  #   }
  # END: FIT DIFFERENT COMBOS
  out.df <- do.call(rbind, out.lst)
  #   print(out.df)
  out.df$bicweight <- BICWeight(out.df$bic)
  
  # get top 3
  out.df <- out.df[order(out.df$bicweight, decreasing = TRUE), ][1:N, ]
  return(out.df)
#   bics.top.i <- order(out.df$bicweight, decreasing = TRUE)[1:3]
#   return(out.df[bics.top.i, ])

#   bics.weight <- BICWeight(out.lst$bics)
#   bics.top.i <- order(bics.weight, decreasing = TRUE)[1:3]
#   fits.top <- [bics.top.i]
#   # get coefficients only: saves memory
#   fits.top <- lapply(fits.top, function(fit) coef(fit))
#   bics.top <- bics.weight[bics.top.i]
#   
#   return(list(fits = fits.top, bics = bics.top))
}

# Source ------------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")

# Load and plot on 12 conditions -------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", verbose = T)

dat.sub <- subset(dat.long, gene == "Nr1d1")
clockgenes <- GetClockGenes()

tiss.test <- c("Liver", "Kidney", "Adr", "BFAT", "Mus")
dat.test <- subset(dat.long, gene %in% clockgenes & tissue %in% c(tiss.test))


# Get design matrices -----------------------------------------------------

experiment <- "experiment"
time <- "time"
exprs <- "exprs"
tissue <- "tissue"
tissues <- as.character(unique(dat.sub$tissue))

form <- GetRhythmicFormula("exprs", "time", c("tissue", "experiment"), with.intercept = FALSE)
des.mat.full <- model.matrix(form, dat.sub)

# tissue.combos.lst <- GetAllCombos(tissues)

# # fit each tissue combo
# fits.all <- lapply(tissue.combos.lst[1:100], function(tiss.sub){
#   if (length(tiss.sub) == 0){
#     rhyth.tiss.str <- NA
#   } else {
#     rhyth.tiss.str <- paste0(tiss.sub, collapse = ",")
#   }
#   des.mat <- SubsetFullDesignMatrix(des.mat.full, rhythmic.tissues = tiss.sub)
#   dat.fits <- dat.test %>%
#     group_by(gene) %>%
#     do(model.fit = lm(exprs ~ 0 + des.mat, .)) %>%
#     mutate(bic = BIC(model.fit), rhyth.tiss = rhyth.tiss.str)
# })
# fits.all

tiss.combos <- GetAllCombos(tiss.test, ignore.full = FALSE)
fits.all <- dat.test %>%
  group_by(gene) %>%
  do(fits = FitCombinations(., tiss.combos))

# size of object
fits.all %>% object.size %>% print(unit = "MB")


# tiss.sub <- combn(tissues, 0)
# des.mat.sub <- SubsetFullDesignMatrix(des.mat.full, tiss.sub[, 1])
# 
# fit <- lm(exprs ~ 0 + des.mats, dat.sub)
# # update through all possibilities, compute BIC
# fit.less <- update(fit, formula. = exprs ~ 0 + des.mat.sub)
# 




