# Jake Yeung
# PCA plotting functions
# Periodogram plotting functions
# Nov 5 2014
# 
PlotLoadings <- function(y, title="Plot title") {
  # Given vector from PCA, plot vector and color by tissue.
  # Maybe give fancy legend
  plot(y, main=title, col=rep(1:12, each=24), type='o') 
}

CalculatePeriodogram <- function(x){
  # Creates periodogram for spectral analysis.
  # 
  # INPUT:
  # x = vector of data you want to check for rhythmicity
  # 
  # OUTPUT:
  # list with $freq being frequencies and periodogram values $periodogram
  #
  # Uses fast fourier transform
  # Adopted from https://onlinecourses.science.psu.edu/stat510/node/71
  
  # BEGIN: do my FFT, extract only relevant ranges
  N <- length(x)
  # Create my "unscaled" periodogram
  FF <- abs(fft(x) / sqrt(N)) ^ 2  # do my FFT. Fast Fourier Transform
  
  # Create my "scaled" periodogram. Scale constant is 4/n
  scale.factor <- 4 / N
  # only need first (N/2) + 1 values of FFT result for periodogram
  P <- scale.factor * FF[1:(N / 2 + 1)]
  P.unscaled <- FF[1:(N / 2 + 1)]
  
  # creates harmonic frequencies from 0 to 0.5 in steps of 1/12 for periodogram.
  # I only need from 0 to 0.5.
  f <- (0:(N / 2) / N)  
  # END: do my FFT, extract only relevant ranges
  
  # Why can't R return two objects? Returning list instead...
  return(list("freq"=f, "p.scaled"=P, "p.unscaled"=P.unscaled))
}

FindMaxFreqs <- function(freq, periodogram, n=5) {
  # Given periodogram and frequency, return frequency at which
  # maximum value of periodogram occurs
  # 
  # Input:
  # f = frequency calculated from CalculatePeriodogram
  # P = periodogram calculated from CalculatePeriodogram
  # n = return the top n frequencies. Default 5
  # 
  # Output:
  # max frequency
  # 
  max.vals <- sort(periodogram, decreasing=TRUE)[1:n]
  max.indices <- match(max.vals, periodogram)
  # Get freqs from indices
  max.freqs <- freq[max.indices]
  
  return(max.freqs)
}

PlotPeriodogram <- function(freq, periodogram, title="Plot title", vline=NA) {
  # Plots periodogram. 
  # 
  # Input: 
  # f = "frequency" calculated from CalculatePeriodogram
  # P = "periodogram" calculated from CalculatePeriodogram
  # 
  # Output:
  # Periodogram plot
  plot(freq, periodogram, type='l', main=title)
}

GetTissueSpecificMatrix <- function(n.tissues, n.timepoints){
  # set my tissue specific factors, matches my response y to a 
  # specific tissue. 12 tissues, 24 time points
  # n.tissues <- 12  # 12 conditions
  # n.timepoints <- 24  # 24 time points
  
  # initialize first tissue vector, then cbind for next tissues
  tissue.factors <- rep(0, n.tissues * n.timepoints)
  tissue.factors[1:n.timepoints] <- rep(1, n.timepoints)
  # print(length(tissue.factors))
  
  # do same for all other tissues, cbind to tissue.factors
  for (c in 2:n.tissues){  # start at 2 because we did 1 already
    start.i <- (c - 1) * n.timepoints + 1  # starts at 25, if c = 2
    end.i <- start.i + n.timepoints - 1  # ends at 48, if c = 2 
    t.fac <- rep(0, n.tissues * n.timepoints)
    t.fac[start.i:end.i] <- rep(1, n.timepoints)
    # print(length(t.fac))
    tissue.factors <- cbind(tissue.factors, t.fac)
  }
  return(tissue.factors)
}

