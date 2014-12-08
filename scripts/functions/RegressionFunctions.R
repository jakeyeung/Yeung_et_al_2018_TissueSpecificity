# Functions involving linear regressions and setting up matrices to get it done.
# used in fit_array_with_exprs.rnaseq.vs.array.R and other similar scripts.
# November 24 2014
# Jake Yeung

InitCoeffMat <- function(row.names, col.names){
  # Prepare empty matrix of NAs with rownames as supplied.
  # Each colname is split into two columns in coeff.mat, a samp1_intercept and samp1_slope
  # 
  # Args
  # row.names: Rownames of coeff.mat
  # colnames: Colnames of coeff.mat
  # 
  # Output
  # coeff.mat: empty matrix of NAs containing rownames and colnames.
  
  coeff.mat <- matrix(NA, nrow=length(row.names), ncol=(length(col.names) * 2))
  # set up row and colnames of coeff.mat
  rownames(coeff.mat) <- row.names
  coeff.mat.colnames <- rep(NA, 2 * length(col.names))
  for (i in 1:length(col.names)){
    col <- col.names[i]
    coeff.mat.colnames[i * 2] <- paste0(col, '_slope')
    coeff.mat.colnames[i * 2 - 1] <- paste0(col, '_intercept')
  }
  colnames(coeff.mat) <- coeff.mat.colnames
  return(coeff.mat)
}

LmGeneTissue <- function(array.subset, rna.seq, row.names, 
                         tissue.names, n.samps){
  # Do lm fit to each tissue, for each gene.
  # 
  # Args:
  # array.subset: exprs of array, has same colnames as rna.seq
  # rna.seq: exprs of rnaseq
  # row.names: gene names
  # tissue.names: Tissue names
  # n.samps: samples in each tissues
  
  coeff.mat <- InitCoeffMat(row.names, tissue.names)
  
  for (gene in row.names){
    for (tissue.i in 1:length(tissue.names)){
      tissue <- tissue.names[tissue.i]
      tissue.i.end <- tissue.i * n.samps
      tissue.i.start <- tissue.i.end - n.samps + 1
      y <- rna.seq[gene, tissue.i.start:tissue.i.end]
      x <- array.subset[gene, tissue.i.start:tissue.i.end]
      fit <- lm(y ~ x)
      intercept <- fit$coefficient[1]
      slope <- fit$coefficient[2]
      coeff.mat[gene, paste0(tissue, '_slope')] <- slope
      coeff.mat[gene, paste0(tissue, '_intercept')] <- intercept
    }
  }
  return(coeff.mat)
}


AdjustArrayToRnaSeq <- function(array.exprs, coeff.mat, tissue.names){
  # init output df
  array.exprs.adjusted <- matrix(NA, nrow=nrow(array.exprs), ncol=ncol(array.exprs),
                                 dimnames=list(rownames(array.exprs),
                                               colnames(array.exprs)))
  for (tissue in tissue.names){
    intercept <- coeff.mat[, paste0(tissue, '_intercept')]
    slope <- coeff.mat[, paste0(tissue, '_slope')]
    tissue.exprs.array <- array.exprs[, grepl(tissue, 
                                              colnames(array.exprs))]
    tissue.exprs.array.normalized <- intercept + slope * tissue.exprs.array
    # write to adjusted exprs
    array.exprs.adjusted[, grepl(tissue, 
                                 colnames(array.exprs.adjusted))] <- tissue.exprs.array.normalized
  }
  return(array.exprs.adjusted)
}

ConstrainedFitWithNoise <- function(gene, array.subset, array.exprs, rna.seq, noise.function,
                                    noise.params){
  # functions
  R.hat <- function(m, a, b) a * (m - b) # RNA to Microarray model.
  S <- function(x) sum((R - R.hat(M, x[1], x[2])) ^ 2 / R.var)  # optimization equation
  
  # fit for every probe a weighted linear regression with weights of 1/var
  R <- rna.seq[gene, ]
  M <- array.subset[gene, ]
  M.full <- array.exprs[gene, ]
  
  a.init <- 1.01  # expect to be > 1
  b.init <- 0.5 * min(M.full)  # background some fraction of min exprs
  
  a.min <- 0
  a.max <- 2^10
  b.min <- 0  # background can't be negative
  b.max <- min(M.full)  # background can't be larger than min exprs
  
  R.var <- noise.function(noise.params, R)
  
  a.b <- optim(c(a.init, b.init), S, method="L-BFGS-B",
               lower=c(a.min, b.min),
               upper=c(a.max, b.max))
  a.hat <- a.b$par[1]
  b.hat <- a.b$par[2]
  conv <- a.b$convergence
  
  R.predict <- R.hat(M.full, a.hat, b.hat)
  
  # add to coeff.mat
  # return(list(a.hat, b.hat, conv))
}

FTestSigmoidLinearModels <- function(fit.lm, fit.complex, pval=1e-10, complex.model="sigmoid"){
  if (!all(is.na(fit.complex))){
    # has sigmoid fit, compare with linear fit with F test
    f.test <- anova(fit.complex, fit.lm)
    f.test.pval <- f.test[["Pr(>F)"]][[2]]
    if (f.test.pval < pval){
      # it is "worth it" to add parameters and fit sigmoid
      myfit <- fit.complex
      fit.used <- complex.model
    } else {
      # just stick with linear
      myfit <- fit.lm
      fit.used <- "lm"
    }
  } else {
    # sigmoid didn't work, suspect the data was not in a shape that allowed
    # good convergence with a sigmoidal function. Stick with linear.
    myfit <- fit.lm
    fit.used <- "lm"
  }
  return(list(myfit=myfit, fit.used=fit.used))
}

FitLmConstraint <- function(M.exprs, R.exprs, weights, 
                            int0, slope0, 
                            intmin, slopemin, 
                            intmax, slopemax){
  # fit model y = b + a * x with constraints.
  # if error, then do y = b + a * x without constraints.
  fit.lm <- tryCatch({
    fit.lm <- nls(M.exprs ~ int + slope * R.exprs,
                  algorithm = "port",
                  start=list(int=int0,
                             slope=slope0),
                  lower=list(int=intmin,
                             slope=slopemin),
                  upper=list(int=intmax,
                             slope=slopemax),
                  weights=weights)
  }, error = function(e){
    warning(e)  # print error message as a warning message.
    fit.lm <- lm(M.exprs ~ R.exprs, 
                 weights=weights)
    return(fit.lm)
  })
  return(fit.lm)
}

PvalFromFit <- function(fit){
  # From lm fit, get a pvalue from beta distribution of R squared value.
  # http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
  # for more details
  # 
  # Args:
  # fit from lm
  # 
  # returns:
  # pval from beta distribution of R squared value.
  p <- summary(fit)$df[1]  # number of parameters in your model
  n.minus.p <- summary(fit)$df[2]  # number of data points minus parameters
  R.squared <- summary(fit)$r.squared[[1]]  # R.squared value of fit
  pval <- pbeta(R.squared, (p - 1)/2, n.minus.p / 2, lower.tail = FALSE, log.p = FALSE)
  return(pval)
}