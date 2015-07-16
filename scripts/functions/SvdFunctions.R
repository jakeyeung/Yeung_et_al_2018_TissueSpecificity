# SvdFunctions.R
# to be used in svd.on.rhythmic.genes.R and other similar scripts
# March 2 2015

# library(plyr)
library(reshape2)
# library(PhaseHSV)
library(ggplot2)

ProjectToFrequency2 <- function(dat, omega, add.tissue=FALSE){
  # simpler than ProjectToFrequency().
  # expect dat to be gene i in tissue c with column names time and exprs
  exprs.transformed <- DoFourier(dat$exprs, dat$time, omega = omega)
  if (add.tissue){    
    tissue <- unique(dat$tissue)
    dat.out <- data.frame(tissue = tissue, exprs.transformed = exprs.transformed)
  } else {
    dat.out <- data.frame(exprs.transformed = exprs.transformed)
  }
  return(dat.out)
}

ProjectToFrequency <- function(dat, my.omega, normalize = TRUE, rhythmic.only = FALSE, method = "ANOVA", pval.cutoff = 5e-3){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # dat: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
  # rhytmic.only: transforms only genes that are rhythmic (by BIC method)
  # 
  # omega in which we are interested.  
  if (rhythmic.only){
    if (IsRhythmic(dat, my.omega, pval.cutoff = pval.cutoff, method = "ANOVA")){
      dat.transformed <- Transform(dat, my.omega, normalize)
    } else {
      dat.transformed <- data.frame(exprs.transformed = 0)
    }
  } else {
    dat.transformed <- Transform(dat, my.omega, normalize)
  }
  return(dat.transformed)
}



Transform <- function(dat, my.omega, normalize = TRUE){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # dat: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
  # 
  # omega in which we are interested.
  
  # if we normalize, get list of omegas.
  # otherwise we just use omega
  if (normalize){
    # t <- sort(unique(dat$time))
    # n.timepoints <- length(t)
    # interval <- t[2] - t[1]
    omegas <- GetOmegas(remove.zero = TRUE)
  } else {
    omegas <- my.omega
  }
  
  transforms <- sapply(omegas, DoFourier, exprs = dat$exprs, time = dat$time)
  my.transformed <- transforms[which(omegas == omega)]  # corresponds to omega of interest
  
  if (normalize){
    # Normalize across omega
    # if median is 0, then set factor to 0, otherwise 1
    factor <- 1
    cutoff <- 5
    jmedian <- median(subset(dat, experiment == "rnaseq")$exprs)
    if (jmedian <= cutoff){
      factor <- 0
    }
    
    my.transformed <- (my.transformed / sqrt(sum(Mod(transforms) ^ 2))) * factor
  } 
  return(data.frame(exprs.transformed = my.transformed))
}

GetInterval <- function(time.vector){
  # Given vector of times (equally spaced time points), return the time interval
  time.sort <- sort(unique(dat$time))
  interval <- time.sort[2] - time.sort[1]
  return(interval)
}

IsRhythmic <- function(dat, my.omega, pval.cutoff = 5e-3, method = "ANOVA"){
  # Test if rhythmic by BIC model selection: fit through all omegas.
  # dat: long format, gene and condition
  # my.omega: omega of interest
  # method = "BIC" or "ftest"
  # pval.cutoff: for method = "ftest"
  
  fit.rhyth <- lm(exprs ~ 0 + experiment + sin(my.omega * time) + cos(my.omega * time), data = dat)
  fit.flat <- lm(exprs ~ 0 + experiment, data = dat)  # intercept only
  
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
                                        data = dat))
    bic.noise.min <- min(sapply(fits.noise, BIC))
    # print(paste(unique(dat$gene), "noise bic:", bic.noise.min, "rhyth bic:", rhyth.bic))
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

# fit.rhyth <- lm(exprs ~ sin(omega * time) + cos(omega * time), data = dat)
# fit.flat <- lm(exprs ~ 1, data = dat)  # intercept only
# # omega.noise <- 2 * pi / (GetInterval(dat$time) * 4)  # period of every two intervals is noise.
# omega.noise <- 2 * pi / 12
# fit.noise <- lm(exprs ~ sin(omega.noise * time) + cos(omega.noise * time), data = dat)
# bic.test <- BIC(fit.rhyth, fit.flat, fit.noise)
# (bic.test)
# # extract rowname of minimum BIC
# chosen.model <- rownames(bic.test)[which(bic.test$BIC == min(bic.test$BIC))]
# if (chosen.model == "fit.rhyth"){
#   dat.transformed <- Transform(dat, my.omega, normalize)
# } else {
#   # if not rhythmic, then do not do transform it.
#   dat.transformed <- data.frame(exprs.transformed = 0)
# }
# } else {
#   dat.transformed <- Transform(dat, my.omega, normalize)
# }

DoFourier <- function(exprs, time, omega, normalize = TRUE){
  p <- exprs %*% exp(1i * omega * time)
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

TemporalToFrequency <- function(dat, period = 24){
  library(plyr)
  omega <- 2 * pi / period
  
  start.time <- Sys.time()
  dat.complex <- lapply(split(dat, dat$tissue), function(x){
    ddply(x, .(gene), ProjectToFrequency2, omega = omega, add.tissue = TRUE)
  }) %>%
    do.call(rbind, .) %>%
    mutate(magnitude = Mod(exprs.transformed)) %>%
    arrange(desc(magnitude))
  print(Sys.time() - start.time)
  
  detach("package:plyr", unload=TRUE)
  library(dplyr)
  return(dat.complex)
}

GetEigens <- function(s.complex, period, comp = 1){
  source("scripts/functions/PlotFunctions.R")
  if (missing(period)){
    period <- 24
  }
  omega <- 2 * pi / period
  
  var.explained <- s.complex$d ^ 2 / sum(s.complex$d ^ 2)
  eigengene <- s.complex$v[, comp]
  eigensamp <- s.complex$u[, comp]
  # rotate to phase of largest magnitude in sample of eigengene
  phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
  rotate.factor <- complex(modulus = 1, argument = phase.reference)
  # rotate eigengene by -phase ref
  eigengene <- eigengene * Conj(rotate.factor)
  # rotate eigensamp by +phase ref
  eigensamp <- eigensamp * Conj(rotate.factor)
  v.plot <- PlotComplex2(eigengene, labels = rownames(s.complex$v), omega = omega, title = paste0("Right singular value ", comp, " (", signif(var.explained[comp], 2), ")"))  
  u.plot <- PlotComplex2(eigensamp, labels = rownames(s.complex$u), omega = omega, title = paste0("Left singular value ", comp, " (", signif(var.explained[comp], 2), ")"))
  return(list(v.plot = v.plot, u.plot = u.plot, eigengene = eigengene, eigensamp = eigensamp))
}