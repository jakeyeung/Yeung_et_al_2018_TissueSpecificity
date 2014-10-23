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

PlotPeriodogram <- function(x, N, title="Plot title", vline=NA) {
  # Creates periodogram for spectral analysis.
  # x is vector of data you want to check for rhythmicity
  # Uses fast fourier transform
  # Adopted from https://onlinecourses.science.psu.edu/stat510/node/71
  
  # BEGIN: do my FFT, extract only relevant ranges
  N <- length(x)
  FF <- abs(fft(x) / sqrt(N)) ^ 2  # do my FFT. Fast Fourier Transform
  # only need first (N/2) + 1 values of FFT result for periodogram
  P <- FF[1:(N / 2 + 1)]
  # creates harmonic frequencies from 0 to 0.5 in steps of 1/12 for periodogram.
  # I only need from 0 to 0.5.
  f <- (0:(N / 2) / N)  
  # END: do my FFT, extract only relevant ranges
  
  # Plot my periodogram
  plot(f, P, type='l', main=title)
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


# Plot PCA ----------------------------------------------------------------

# Loop to plot scatter plot of PCA i versus PCA i + 1
for (x_comp in 1:5) {

  y_comp <- x_comp + 1  # plot x PCA component i on x-axis, PCA component i+1 on y-axis
  

  # Show PCA plot of tissues ----------------------------------------------------
  
  # TWO WAYS OF COLORING SAMPLES IN PCA
  
  # 1.
  #   # Color by time point
  #   textplot(dat_pca$x[, x_comp], dat_pca$x[, y_comp], newnames, cex=0.7, col=rep(1:24, 12), 
  #           main=paste("Component", x_comp, "vs.", y_comp))
  
  # 2. 
  # Color by tissue
  textplot(dat_pca$x[, x_comp], dat_pca$x[, y_comp], newnames, cex=0.7, col=rep(1:12, each=24), 
           main=paste("Component", x_comp, "vs.", y_comp))
}

# Which vector loadings have oscillating components? ----------------------

N <- nrow(dat_pca$x)  # number of samples.
T <- 24  # 24 hours in a period

# Fit PCA component of interest: loop to try many different PCAs
# user changeable range
for (pca_vector in 8:20) {
  # Create response vector, which is loadings
  
  y <- dat_pca$x[, pca_vector]
  
  # BEGIN: plot periodograms to see which frequency has high activity
  PlotPeriodogram(y, title=paste("Periodogram for PCA component:", pca_vector))
  # spec.pgram(y, log="no", main=pca_component_string)
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
  
  # Print some statements to describe fit
  print("*********************************")
  print(paste0(pca_vector))
  cat("\n")
  print(anova(fit))
  cat("\n")
}
 


