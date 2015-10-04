BICFromLmFit <- function(coefficients, residuals){
  # from vector of coefs and residuals, calculate BIC
  n <- length(residuals)
  RSS <- sum(residuals ^ 2)
  k <- length(coefficients)
  criterion <- n * log(RSS / n) + k * log(n)
  return(criterion)
}

FitModels <- function(dat.gene, my_mat, get.criterion = "BIC", normalize.weights = TRUE){
  # Fit many models with lm.fit() which is faster than lm()
  weight.sum <<- 0  # track
  fits <- lapply(my_mat, function(mat, my_mat, model.selection = "BIC"){
    fit <- lm.fit(y = dat.gene$exprs, x = mat)
    if (model.selection == "BIC"){
      criterion <- BICFromLmFit(fit$coefficients, fit$residuals)
#       n <- length(fit$residuals)
#       RSS <- sum(fit$residuals ^ 2)
#       k <- length(fit$coefficients)
#       criterion <- -2 * log(RSS / n) + k * log(n)
      weight <- exp(-0.5 * criterion)
    } else {
      weight <- NA
    }
    weight.sum <<- weight.sum + weight
    return(list(fit = fit$coefficients, residuals = fit$residuals, weight = weight))
  }, my_mat)
  if (normalize.weights){
    fits <- lapply(fits, function(fit){
      fit$weight.norm <- fit$weight / weight.sum
      return(fit)
    })
  }
  return(fits)
}

GetSelectionCriterion <- function(fits, model.selection = "BIC"){
  # given list of fits, calculate BIC weights
  # get coefficients and residuals, get model selection by BIC or AIC
  if (model.selection == "BIC"){
    # https://en.wikipedia.org/wiki/Bayesian_information_criterion#Definition
    # BIC = -2 log(RSS/n) + k * log(n)
    bics <- unlist(lapply(fits, function(fit){
      n <- length(fit$residuals)
      RSS <- sum(fit$residuals ^ 2)
      k <- length(fit$coefficients)
      criterion <- -2 * log(RSS / n) + k * log(n)
      bicweight <- exp(-0.5 * criterion)
      }))
    # normalize
    bics.weight <- bics / sum(bics)
  } else {
    warning("Only BIC is currently implemented.")
  }
}

MakeRhythmicDesignMatrices <- function(dat.gene, w = 2 * pi / 24, simplify=FALSE){
  # dat.gene: long format of gene expression, time and
  # conditions. Can include additional factors such as 
  # experiment.
  # simplify: return only design matrix rather than full matrix containing meta data (FALSE is debug mode)
  
  tissues <- unique(as.character(dat.gene$tissue))
  
  tiss.combos <- GetAllCombos(tissues, ignore.full = FALSE)
  my_mat.queue <- new.queue()
  
  # BEGIN: init with flat model
  des.mat.flat <- GetFlatModel(dat.gene)
  # get rhythmic parameters which will be used for adding later: hash structure has fast lookup
  des.mat.sinhash <- GetSinCombos(dat.gene, w, tissues, tiss.combos)
  des.mat.coshash <- GetCosCombos(dat.gene, w, tissues, tiss.combos)
  
  rhyth.tiss <- list(character(0))  # needs to track shared and independent parameters, e.g.: c("Liver,Kidney", "Adr") no duplicates allowed
  # n.rhyth <- NRhythmicFromString(rhyth.tiss)  # number of independent rhythmic parameters perhaps? 
  n.rhyth <- NRhythmicFromVector(rhyth.tiss)  # do that later? naw faster if we do it now
  complement <- FilterCombos(tiss.combos, rhyth.tiss)
  des.mat.list <- list(mat=des.mat.flat, rhyth.tiss=rhyth.tiss, n.rhyth=n.rhyth, complement = complement)
  # END: init with flat model
  
  # load up my queue
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
  print(paste("Number of matrices generated:", n.mat.submitted))
  des.mats.list <- des.mats$as.list()
  if (simplify){
    des.mats.list <- lapply(des.mats.list, function(x) return(x$mat))
  }
  return(des.mats.list)
}

MakeModelName <- function(rhyth.tiss.lst, delim = ";"){
  # From a list of rhythmic tissues, generate a hash key to track models
  # therefore: Adr,Kidney;Liver means Adr,Kidney same param, Liver independent param
  lst <- sort(unlist(rhyth.tiss.lst))
  return(paste(lst, collapse = delim))
}

AddRhythmicColumns <- function(des.mat.sinhash, des.mat.coshash, tiss.key){
  # cbind rhythnic columns, rename to track which tissues share this parameter
  col.new <- cbind(des.mat.sinhash[[tiss.key]], des.mat.coshash[[tiss.key]])
  col.new.namebase <- paste(tiss.key, sep = ",")
  colnames(col.new) <- c(paste0(col.new.namebase, ":sin(w * time)"), paste0(col.new.namebase, ":cos(w * time)"))
  return(col.new)
}

GetSinCombos <- function(dat.gene, w, tissues, combos){
  des.mat.sin <- GetRhythModel.sin(dat.gene, w)
  des.mat.sin.hash <- MakeHash(des.mat.sin, tissues, combos)
  return(des.mat.sin.hash)
}

GetCosCombos <- function(dat.gene, w, tissues, combos){
  des.mat.cos <- GetRhythModel.cos(dat.gene, w)
  des.mat.cos.hash <- MakeHash(des.mat.cos, tissues, combos)
  return(des.mat.cos.hash)
}

NRhythmicFromString <- function(rhyth.tiss, jsep = ","){
  # How many tissues are rhythmic given a comma separated stringS
  # "" suggests 0
  if (length(rhyth.tiss) == 0) return(0)
  if (rhyth.tiss == ""){
    return(0) 
  } else {
    return(length(strsplit(rhyth.tiss, jsep)[[1]]))
  }
}

NRhythmicFromVector <- function(rhyth.tiss.vec, jsep = ","){
  # if c("Adr,Liver", "Kidney") output 3 rhythmic tisues
  n.rhyth <- 0
  for (tiss in rhyth.tiss.vec){
    n.rhyth <- n.rhyth + NRhythmicFromString(tiss)
  }
  return(n.rhyth)
}

FilterCombos <- function(combos, combos.sub){
  # Filter out combos given a combos.sub, filters
  # any combo that contains even one of the elements in combo.sub
  # always removes an empty set no matter what
  combos.subl <- sapply(combos, function(comb){
    if (length(comb) == 0) return(FALSE)  # empty set
    inrhyth <- TRUE  # say it is true, loop through conditions to set to false if it is a duplicate
    for (ci in comb){
      if (ci %in% combos.sub){
        inrhyth <- inrhyth * FALSE
      } 
    }
    return(as.logical(inrhyth))
  })
  return(combos[combos.subl])
}

GetRhythmicFormula <- function(exprs, time, intercepts, w = 2 * pi / 24, with.intercept = FALSE){
  # Get formula for full rhythmic model
  # tissue: genotype, or tissues
  # x: time
  # y: expression
  # experiment: column name of experiment (can be blank I think "")
  # with.intercept: if you want tissue and experiments to be your intercepts, keep it False
  
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

GetFlatModel <- function(dat.gene){
  des.mat.flat <- model.matrix(exprs ~ 0 + tissue + experiment, dat.gene)
}

GetRhythModel <- function(dat.gene, w = 2 * pi / 24){
  des.mat.rhyth <- model.matrix(exprs ~ 0 + tissue:sin(w * time) + tissue:cos(w * time), dat.gene)
}

MakeHash <- function(des.mat.rhyth, tissues, combos){
  # return as hash object requires hash library
  tissues <- unique(as.character(dat.gene$tissue))
  # init with single tissues
  des.mat.hash <- hash()
  for (tiss.i in seq(tissues)){
    tiss <- tissues[tiss.i]
    des.mat.hash[[tiss]] <- des.mat.rhyth[, colnames(des.mat.rhyth)[grepl(tiss, colnames(des.mat.rhyth))]]
  }
  # add combos only if they are composed of >1 tissues
  for (comb in combos){
    if (length(comb) <= 1) next
    # init
    vec <- numeric(length = nrow(des.mat.rhyth))
    for (tiss in comb){
      vec <- vec + des.mat.hash[[tiss]]
    }
    # put combo into hash as comma separated string
    comb.key <- paste(comb, collapse = ",")
    des.mat.hash[[comb.key]] <- vec
  }
  return(des.mat.hash)
}


GetRhythModel.sin <- function(dat.gene, w = 2 * pi / 24){
  des.mat.rhyth <- model.matrix(exprs ~ 0 + tissue:sin(w * time), dat.gene)
  return(des.mat.rhyth)
}

GetRhythModel.cos <- function(dat.gene, w = 2 * pi / 24){
  des.mat.rhyth <- model.matrix(exprs ~ 0 + tissue:cos(w * time), dat.gene)
}

GetRhythmicFormula.Shared <- function(exprs, time, rhythmic.factor, intercepts, w = 2 * pi / 24, with.intercept = FALSE){
  # Get formula for full rhythmic model
  # exprs: expression col
  # time: time
  # rhythmic.factor: factor describing if rhythmic or not. 
  # with.intercept: if you want tissue and experiments to be your intercepts, keep it False
  
  # get rhythmic terms
  sinterm <- paste0(rhythmic.factor, " : ", "sin(", w, " * ", time, ")")
  costerm <- paste0(rhythmic.factor, " : ", "cos(", w, " * ", time, ")")
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
    if (is.null(ncol(tissues.combo))){
      # only one choice, put that choice into list
      tissue.combos.lst[[model.i]] <- tissues.combo
      model.i <- model.i + 1
    }
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
  # fit.full <- lm(exprs ~ 0 + des.mat.sub, data = dat.gene)
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
  # gc()  # too slow
  return(out.df)
}

FitRhythmicDesMat <- function(dat, des.mat){
  # Fit rhythmic with design matix
  mod <- lm(exprs ~ 0 + des.mat, data = dat)
  bic <- BIC(mod)
  return(list(fit = coef(mod), bic = bic))
}