# Jake Yeung
# explore_data.R
# initial exploration of data in microarray across tissues

library(ggplot2)
library(reshape)
library(wordcloud)  # for showing text without jumbling


# Define my functions -----------------------------------------------------

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

ShortenSampNames <- function(long.names, show="tissue") {
  # Function to shorten sample names of a certain form. 
  # User options to determine how much to shorten
  #
  # INPUT: 
  # [GSM1321058_BS58_MoGene1.0ST.CEL, GSM1321058_BS59_MoGene1.0ST.CEL ..., ]
  # 
  # OUTPUT:
  # [BS58, BS59, ..., ]  if show == "tissue.time"
  # or
  # [BS, BS, ..., ]  if show == "tissue"
  # or
  # [58, 59, ..., ] if show == "time"
  # 
  # Whether you want the number afterwards is user defined.
  
  long.names <- strsplit(long.names, '_')
  short.names <- sapply(long.names, function(x) {
    sample.id <- x[2]    # sample name with time point
    
    # THREE WAYS OF SHOWING SAMPLES NAMES:
    
    if (show == "tissue.time") {
      # 1.
      # show sample name with timepoint e.g. Liver34
      samp.name <- sample.id 
    }
    else if (show == "tissue") {
      # 2.
      # show opnly tissue component in sample name e.g. "Adr", "Liver"
      samp.name <- substring(sample.id, 1, nchar(sample.id) - 2)      
    }
    else if (show == "time") {
      # 3.
      # show only time component in sample name e.g. 34
      samp.name <- substring(sample.id, nchar(sample.id) - 1, nchar(sample.id)) 
    }
  })
  return(short.names)
}

MeanCenterAcrossGroups <- function(x, n.per.group=24) {
  # Input:
  # x = vector x of inputs across all samples, contains n.per.group
  # elements per group. 
  # Vector is ordered such that all elements
  # in group i are together. 
  # Counting n.per.group would
  # reveal each group. 
  # 
  # Output:
  # x.mean.center <- calculate mean of each group (i.e. mean
  # of every group of n.per.group) and subtract each element
  # in that group by its mean.
  
  N <- length(x)
  x.mean.center <- rep(NA, length(x))
  # BEGIN SANITY CHECK
  # check that N / n.per.group gives no remainder
  if (N %% n.per.group != 0) {
    warning("Uneven number of elements per group.")
  }
  n.groups <- N / n.per.group  # should be an integer
  # END SANITY CHECK
  

  for (i in 0:n.groups - 1) {
    # loop from 0 to n.groups - 1 so that our start index
    # loops like 1, 1 + n.per.group, 1 + 2 * n.per.group
    # for example start.index = 1, 25, 49 ...
    start.index <- i * n.per.group + 1
    end.index <- start.index + n.per.group - 1
    # get group e.g. x[1:24]
    group <- x[start.index:end.index]
    # subtract each element in group by mean(group)
    group.mean.center <- group - mean(group)
    # put into x.mean.center
    x.mean.center[start.index:end.index] <- group.mean.center
  }
  
  return(x.mean.center)
}


# Load data, log transform ------------------------------------------------

# define dirs
data_dir <- "microarray_data"
fname <- "hogenesch_2014_rma.txt"    # data reprocessed by RMA package

# load data
data_path <- file.path(data_dir, fname)
print(paste("Reading data from,", data_path, "May take a few a minutes."))
dat <- read.table(data_path)
print("Read data to memory.")

# log2 transform
dat_original <- dat
dat <- log2(dat)

# make more meaningful sample names
# user changeable paramters: show="tissue.time" | "tissue" | "time"
colnames(dat) <- ShortenSampNames(colnames(dat), show="tissue.time")


# Plot histogram and quantiles ----------------------------------------------------------

# Optional: can comment out, takes time and memory to "melt" data
# dat_melt <- melt(dat)
# 
# ggplot(data=dat_melt, aes(x=value)) +
#   geom_density()
# 
# p <- ggplot(dat_melt, aes(factor(variable), value))
# p + geom_boxplot()


# Calculate PCA and Screeplot ---------------------------------------------

dat_standardize <- as.data.frame(scale(dat))

dat_pca <- prcomp(t(dat_standardize))

screeplot(dat_pca, type="lines", npcs = min(287, length(dat_pca$sdev)), log="y")


# Plot PCA: tissue components ----------------------------------------------------

# Loop to plot scatter plot of PCA i versus PCA i + 1
for (x_comp in 1:5) {

  y_comp <- x_comp + 1  # plot x PCA component i on x-axis, PCA component i+1 on y-axis
  

  # Show PCA plot of tissues ----------------------------------------------------
  
  # TWO WAYS OF COLORING SAMPLES IN PCA: by time or by tissue
  
  
  # 1.
  #   # Color by time point
  #   colors.by.time <- rep(1:24, 12)
  #   textplot(dat_pca$x[, x_comp], dat_pca$x[, y_comp], newnames, cex=0.7, col=colors.by.time, 
  #           main=paste("Component", x_comp, "vs.", y_comp))
  
  # 2. 
  # Color by tissue
  colors.by.tissue <- rep(1:12, each=24)
  textplot(dat_pca$x[, x_comp], dat_pca$x[, y_comp], newnames, cex=0.7, col=colors.by.tissue, 
           main=paste("Component", x_comp, "vs.", y_comp))
  
}


# Plot PCA: circadian components ---------------------------------------------


# Loop to plot scatter plot of PCA i versus PCA i + 1
# FELIX: this is for the 'circadian; components
for (x_comp in 14:19) {
  
  y_comp <- x_comp + 1  # plot x PCA component i on x-axis, PCA component i+1 on y-axis
  
  
  # Show PCA plot of tissues ----------------------------------------------------
  
  # TWO WAYS OF COLORING SAMPLES IN PCA
  
  # 1.
  #   # Color by time point
  colors=hsv(rep(c(1:12)/12, 24), 1, 1)

  ii=1:nrow(dat_pca$x)
#   ii=which(newnames=="Liver")

  newnames.2=paste(newnames, rep(c(0:23)*2, 12))
  
  r=max(range(c(dat_pca$x[ii, x_comp], dat_pca$x[ii, y_comp]))) * 1.2

  textplot(dat_pca$x[ii, x_comp], dat_pca$x[ii, y_comp], newnames.2[ii], cex=0.7, col=colors, 
           main=paste("Component", x_comp, "vs.", y_comp), xlim=c(-r,r), ylim=c(-r,r))
}

# Which vector loadings have oscillating components? ----------------------

# Optional: can MeanCenter the vector y to remove cross-tissue differences
# and focus only on oscillations.

N <- nrow(dat_pca$x)  # number of samples.
T <- 24  # 24 hours in a period

# Fit PCA component of interest: loop to try many different PCAs
# user changeable range
for (pca_vector in 14:20) {
  # Create response vector, which is loadings
  
  y <- dat_pca$x[, pca_vector]
  
  # Optional: 
  # y <- MeanCenterAcrossGroups(y)
  
  # BEGIN: plot periodograms to see which frequency has high activity
  freq.and.periodogram <- CalculatePeriodogram(y)  # returns a list
  freq <- freq.and.periodogram$freq
  periodogram <- freq.and.periodogram$p.scaled
  periodogram.unscaled <- freq.and.periodogram$p.unscaled
  
  # Calculate top 5 frequencies
  max.freqs <- FindMaxFreqs(freq, periodogram)
  
  PlotPeriodogram(freq, periodogram, title=paste("Periodogram for PCA component:", pca_vector))
  # add vertical line at max frequency
  max.f <- max.freqs[1]
  # calculate period from frequency
  max.T <- (1 / max.f) * 2  # multiply by 2 because samples are every 2 hours 
  abline(v=max.f, col='blue', lwd=2)
  # add text to show freq and period.
  # x offset of +0.02 so you can see the text
  text(max.f + 0.02, 0, paste0("T=", signif(max.T, digits=3), "hrs"))
  
  PlotLoadings(y, title=paste("Vector Loadings for PCA component:", pca_vector))
  # END: plot periodograms to see which frequency has high activity
  
  # BEGIN: Linear fit for period of 24 hours
  # Create sequence of t = [0, 2, 4, ... (N - 1)]
  t <- seq(0, 2 * N - 1, 2)
  
  # set my angular frequency
  omega <- (2 * pi) / (T)
  
  # fit my lm using cos and sin with angular frequency
  fit <- lm(y ~ cos(omega * t) + sin(omega * t))
  # END: Linear fit for period of 24 hours
  
  # Print some statementsy to describe fit
  cat("*********************************\n")
  cat(paste0("PCA ", pca_vector, "\n"))
  cat("\n")
  
  cat("Max freqs\n")
  print(max.freqs)
  
  print(anova(fit))
  cat("\n")
}


# Limma: multiple intercept fit -------------------------------------------

# fit linear model: with multiple intercept but shared oscillating factor

# Optional: can MeanCenter the vector y to remove cross-tissue differences
# and focus only on oscillations.

N <- nrow(dat_pca$x)  # number of samples.
T <- 24  # 24 hours in a period

# Fit PCA component of interest: loop to try many different PCAs
# user changeable range
for (pca_vector in 1:20) {
  # Create response vector, which is loadings
  
  y <- dat_pca$x[, pca_vector]
  
  # Optional: 
  # y <- MeanCenterAcrossGroups(y)
  
  # BEGIN: plot periodograms to see which frequency has high activity
  freq.and.periodogram <- CalculatePeriodogram(y)  # returns a list
  freq <- freq.and.periodogram$freq
  periodogram <- freq.and.periodogram$p.scaled
  periodogram.unscaled <- freq.and.periodogram$p.unscaled
  
  # Calculate top 5 frequencies
  max.freqs <- FindMaxFreqs(freq, periodogram)
  
  PlotPeriodogram(freq, periodogram, title=paste("Periodogram for PCA component:", pca_vector))
  # add vertical line at max frequency
  max.f <- max.freqs[1]
  # calculate period from frequency
  max.T <- (1 / max.f) * 2  # multiply by 2 because samples are every 2 hours 
  abline(v=max.f, col='blue', lwd=2)
  # add text to show freq and period.
  # x offset of +0.02 so you can see the text
  text(max.f + 0.02, 0, paste0("T=", signif(max.T, digits=3), "hrs"))
  
  PlotLoadings(y, title=paste("Vector Loadings for PCA component:", pca_vector))
  # END: plot periodograms to see which frequency has high activity
  
  # BEGIN: Linear fit for period of 24 hours
  # Create sequence of t = [0, 2, 4, ... (N - 1)]
  t <- seq(0, 2 * N - 1, 2)
  
  # set my angular frequency
  omega <- (2 * pi) / (T)
  
  # set my tissue specific factors, matches my response y to a 
  # specific tissue. 12 tissues, 24 time points
  n.tissues <- 12  # 12 conditions
  n.timepoints <- 24  # 24 time points
  
  # initialize first tissue vector, then cbind for next tissues
  tissue.factors <- rep(0, n.tissues * n.timepoints)
  tissue.factors[1:24] <- rep(1, n.timepoints)
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
  
  print(dim(tissue.factors))
  # fit my lm using cos and sin with angular frequency
  fit <- lm(y ~ 0 + tissue.factors + cos(omega * t) + sin(omega * t))
  # END: Linear fit for period of 24 hours
  
  # Print some statementsy to describe fit
  cat("*********************************\n")
  cat(paste0("PCA ", pca_vector, "\n"))
  cat("\n")
  
  print(fit$coefficients)
  print(anova(fit))
  cat("\n")
}


 


