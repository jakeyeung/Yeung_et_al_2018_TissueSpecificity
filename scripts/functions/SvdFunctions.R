# SvdFunctions.R
# to be used in svd.on.rhythmic.genes.R and other similar scripts
# March 2 2015

Transform <- function(df, my.omega, normalize = TRUE){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # df: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
  # 
  # omega in which we are interested.
  
  # if we normalize, get list of omegas.
  # otherwise we just use omega
  if (normalize){
    # t <- sort(unique(df$time))
    # n.timepoints <- length(t)
    # interval <- t[2] - t[1]
    omegas <- GetOmegas(remove.zero = TRUE)
  } else {
    omegas <- my.omega
  }
  
  transforms <- sapply(omegas, DoFourier, exprs = df$exprs, time = df$time)
  my.transformed <- transforms[which(omegas == omega)]  # corresponds to omega of interest
  
  if (normalize){
    # Normalize across omega
    # if median is 0, then set factor to 0, otherwise 1
    factor <- 1
    cutoff <- 5
    jmedian <- median(subset(df, experiment == "rnaseq")$exprs)
    if (jmedian <= cutoff){
      factor <- 0
    }
    
    my.transformed <- (my.transformed / sqrt(sum(Mod(transforms) ^ 2))) * factor
  } 
  return(data.frame(exprs.transformed = my.transformed))
}

GetInterval <- function(time.vector){
  # Given vector of times (equally spaced time points), return the time interval
  time.sort <- sort(unique(df$time))
  interval <- time.sort[2] - time.sort[1]
  return(interval)
}

IsRhythmic <- function(df, my.omega, pval.cutoff = 5e-3, method = "ANOVA"){
  # Test if rhythmic by BIC model selection: fit through all omegas.
  # df: long format, gene and condition
  # my.omega: omega of interest
  # method = "BIC" or "ftest"
  # pval.cutoff: for method = "ftest"
  
  fit.rhyth <- lm(exprs ~ 0 + experiment + sin(my.omega * time) + cos(my.omega * time), data = df)
  fit.flat <- lm(exprs ~ 0 + experiment, data = df)  # intercept only
  
  if (method == "BIC"){
    omegas.all <- GetOmegas(remove.zero = TRUE)
    # remove also my.omega to get omegas representing "noise"
    omegas.noise <- omegas.all[which(omegas.all != my.omega)]
    
    bic.test <- BIC(fit.rhyth, fit.flat)
    
    chosen.model <- rownames(bic.test)[which(bic.test$BIC == min(bic.test$BIC))]
    
    if (chosen.model == "fit.flat"){
      # if flat, no need to check for noise.
      return(FALSE)
    }
    # chosen model is fit.rhyth, but is it a noisy gene?
    # Check if noise components have an even better fit than fit.rhyth
    rhyth.bic <- bic.test["fit.rhyth", "BIC"]
    # fit noise
    fits.noise <- lapply(omegas.noise, 
                         function(w) lm(exprs ~ sin(w * time) + cos(w * time), 
                                        data = df))
    bic.noise.min <- min(sapply(fits.noise, BIC))
    # print(paste(unique(df$gene), "noise bic:", bic.noise.min, "rhyth bic:", rhyth.bic))
    if (rhyth.bic < bic.noise.min){
      return(TRUE)
    } else {
      return(FALSE)
    } 
  } else if (method == "ANOVA"){
    f.test <- anova(fit.flat, fit.rhyth)
    pval <- f.test[["Pr(>F)"]][[2]]
    if (is.nan(pval)){
      # if gene is all flat, then pval is NaN, force it to be 1 in this case.
      pval <- 1
    }
    if (pval < pval.cutoff){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

ProjectToFrequency <- function(df, my.omega, normalize = TRUE, rhythmic.only = FALSE, method = "ANOVA", pval.cutoff = 5e-3){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # df: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
  # rhytmic.only: transforms only genes that are rhythmic (by BIC method)
  # 
  # omega in which we are interested.
  
  if (rhythmic.only){
    if (IsRhythmic(df, my.omega, pval.cutoff = pval.cutoff, method = "ANOVA")){
      df.transformed <- Transform(df, my.omega, normalize)
    } else {
      df.transformed <- data.frame(exprs.transformed = 0)
    }
  } else {
    df.transformed <- Transform(df, my.omega, normalize)
  }
  return(df.transformed)
}

# fit.rhyth <- lm(exprs ~ sin(omega * time) + cos(omega * time), data = df)
# fit.flat <- lm(exprs ~ 1, data = df)  # intercept only
# # omega.noise <- 2 * pi / (GetInterval(df$time) * 4)  # period of every two intervals is noise.
# omega.noise <- 2 * pi / 12
# fit.noise <- lm(exprs ~ sin(omega.noise * time) + cos(omega.noise * time), data = df)
# bic.test <- BIC(fit.rhyth, fit.flat, fit.noise)
# (bic.test)
# # extract rowname of minimum BIC
# chosen.model <- rownames(bic.test)[which(bic.test$BIC == min(bic.test$BIC))]
# if (chosen.model == "fit.rhyth"){
#   df.transformed <- Transform(df, my.omega, normalize)
# } else {
#   # if not rhythmic, then do not do transform it.
#   df.transformed <- data.frame(exprs.transformed = 0)
# }
# } else {
#   df.transformed <- Transform(df, my.omega, normalize)
# }

DoFourier <- function(exprs, time, omega, normalize = TRUE){
  p <- exprs %*% exp(-1i * omega * time)
  if (normalize){
    p <- p / length(time)  # normalize by number of datapoints 
  }
  return(p)
}

GetOmegas <- function(n.timepoints=24, interval=2, remove.zero=FALSE){
  # Get omegas from a time series
  # 
  # remove.zero: removes first element (omega = 0).
  # Useful for normalizing Fourier when I don't want omega = 0 case
  
  # Output:
  # list of omegas: can be used for fourier transform
  
  t <- seq(from = 0, by = 1, length.out = n.timepoints)  # will account for interval later
  
  # get harmonic frequencies
  # this creates harmonic frequencies from 0 to .5 in steps of 1/n.timepoints
  
  mid.point <- length(t) / 2 + 1
  freqs <- t[1:mid.point] / n.timepoints
  
  # Convert to period before adjusting for interval
  T.unadj <- 1 / freqs
  
  # Account for interval
  T <- T.unadj * interval
  
  # Calculate omegas
  omegas <- 2 * pi / T
  
  if (remove.zero){
    omegas <- omegas[2:length(omegas)]
  }
  return(omegas)
}

ConvertLongToWide <- function(long.df, measurement.var = "exprs.transformed"){
  wide.df <- dcast(dat.proj, gene ~ tissue, value.var = measurement.var)
  # first row is gene name, let's make them rowname and remove first column.
  
  rownames(wide.df) <- wide.df$gene
  
  wide.df <- subset(wide.df, select = -c(gene))
  
  return(wide.df)
}

PlotComponentOriginal <- function(orig.dat, s, jgene, component, show.original = TRUE){
  component.mat <- s$d[component] * OuterComplex(s$u[, component, drop = FALSE], t(s$v[, component, drop = FALSE]))
  if (show.original){
    PlotComplex(as.matrix(orig.dat[jgene, ]), labels = colnames(orig.dat), main = paste(jgene, "original"), add.text.plot = FALSE) 
  }
  PlotComplex(as.matrix(component.mat[jgene, ]), labels = colnames(component.mat), main = paste(jgene, "component", component), add.text.plot = FALSE)
}

ProjectOnVector <- function(input.complex, basis.complex){
  # Given an input complex number, project onto a basis complex number using
  # simple trigonometry or dot product
  # 
  # Returns length of vector after projection. >0 indicates along basis vector.
  # <0 indicates anti-correlation with basis vector.
  
  if (is.data.frame(input.complex)){
    input.complex <- as.matrix(input.complex)
  }
  if (is.data.frame(basis.complex)){
    basis.complex <- as.matrix(basis.complex)
  }
  phase.diff <- Arg(input.complex) - Arg(basis.complex)
  projected.mod <- Mod(input.complex) * cos(phase.diff)
  return(projected.mod)
}