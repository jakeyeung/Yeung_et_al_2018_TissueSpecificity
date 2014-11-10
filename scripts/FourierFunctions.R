# Jake Yeung
# Nov 10 2014
# functions related to Mr. Fourier

CalculatePeriodogram <- function(x, is.matrix=FALSE){
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
  if (is.matrix==FALSE){
    N <- length(x)  
  }
  else{
    N <- ncol(x)
  }
  # Create my "unscaled" periodogram
  
  if (is.matrix==FALSE){
    FF <- abs(fft(x) / sqrt(N)) ^ 2  # do my FFT. Fast Fourier Transform
  }
  else {
    FF <- t(mvfft(t(x)))
    # FF <- t(abs(mvfft(t(x))))
    # FF <- t(abs(mvfft(t(x)) / sqrt(N)) ^ 2)  # do my FFT. Fast Fourier Transform
  }
  # Create my "scaled" periodogram. Scale constant is 1/n
  scale.factor <- 1 / N
  # only need first (N/2) + 1 values of FFT result for periodogram
  
  if (is.matrix==FALSE){
    P <- scale.factor * FF[1:(N / 2 + 1)]
    P.unscaled <- FF[1:(N / 2 + 1)]    
  }
  else {
    P <- scale.factor * FF[, 1:(N / 2 )]
    P.unscaled <- FF[, 1:(N / 2 )]  
  }
  
  
  # creates harmonic frequencies from 0 to 0.5 in steps of 1/N for periodogram.
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

ProjectToPeriodicTime <- function(Y , N.TISSUES, N.TIMEPTS, INTERVAL, OMEGA, col.names){
  # Project matrix of times to frequency domain
  # ARGS:
  #   Y: matrix of expression. Rows are genes. Columns contain tissues and time.
  #      In this case, expect tissues clustered together, ordered by time.
  #   col.names: column names of output Y's projected onto time. FALSE means no col.names
  #   N.TISSUES: number of tissues
  #   N.TIMEPTS: number of time points per tissue
  #   INTERVAL: interval between time points. e.g. 2 if sampled every 2hrs
  #   OMEGA: 2 * pi / PERIOD. Angular frequency. If omega = 0, matrix Y is 
  #   normalized by time components (not oscillating)
  # 
  # RETURNS:
  #   Y.time.projected: matrix of expression, projected onto the temporal axis.
  
  # track number of genes for dimension purposes
  N.GENES <- nrow(Y)
  # get times vector
  times.vec <- seq(length.out = N.TIMEPTS, by = INTERVAL)
  
  # init output matrix
  Y.time.projected <- matrix(NA, nrow=N.GENES, ncol=N.TISSUES)
  
  # identify row and colnames
  rownames(Y.time.projected) <- rownames(Y)  # same row names
  colnames(Y.time.projected) <- col.names  # from args  
  
  # BEGIN: project onto temporal axis
  for (i in 1:N.TISSUES){
    # get tissue.i across time
    index.start <- (i - 1) * N.TIMEPTS + 1  # 1, 25, 49...
    index.end <- i * N.TIMEPTS
    Y.tissue.i <- Y[, index.start:index.end]  # all genes
    Y.fft <- CalculatePeriodogram(Y.tissue.i, is.matrix=TRUE)
    
    # WHICH FREQUENCY TO EXTRACT? EXTRACT FREQUENCY AT OMEGA.
    freq <- OMEGA / (2 * pi)
    
    # CalculatePeriodogram frequencies are in steps from
    # (0:n.steps/2) by steps 1 / n.steps 
    # we need to divide frequency by INTERVAL to account for sampling interval
    freq.adj <- Y.fft$freq / INTERVAL
    freq.i <- which(freq.adj == freq)
    
    Y.time.projected[, i] <- as.matrix(Y.fft$p.scaled)[, freq.i]
    # print(i)
    # print(Y.tissue.i)
    # project tissues onto temporal axis
    # if OMEGA = 0, it is equivalent to getting average of Y.tissue.i
    # T <- exp(-1i * OMEGA * times.vec) / length(times.vec)
    # Y.time.projected[, i] <- Y.tissue.i %*% T
    # print(T)
    # print(Y.time.projected[, i])
  }
  return(Y.time.projected)
}

# Test function works -----------------------------------------------------

# test ProjectToPeriodicTIme
#
ROWS <- 3
COLS <- 288
n.timepts <- 24  # 24 time points per tissue
n.tiss <- COLS / n.timepts
P <- 24  # 24 hour period
w <- 2 * pi / P
# w <- 0
# use a cosine or sine function with period of 24 hours
t <- seq(from=0, by=2, length.out=COLS)  # 0 to 48 hours sampled every 2 hrs. 12 tissues
y <- 3 * sin(w * t)
# y <- rep(1:n.tiss, each=n.timepts)
Y <- matrix(y, nrow=ROWS, ncol=COLS, byrow = TRUE)
out.colnames <- make.unique(rep('COL', n.tiss))
rownames(Y) <- make.unique(rep('ROW', ROWS))
(Y.t <- ProjectToPeriodicTime(Y, N.TISSUES=n.tiss, N.TIMEPTS=n.timepts, INTERVAL=2, OMEGA=w, out.colnames))
(Mod(Y.t))
