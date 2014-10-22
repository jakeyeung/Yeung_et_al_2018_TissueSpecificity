# Jake Yeung
# explore_data.R
# initial exploration of data in microarray across tissues

library(ggplot2)
library(reshape)
library(wordcloud)  # for showing text without jumbling


# Define my functions -----------------------------------------------------

plot_loadings <- function(y, title="Plot title"){
  # Given vector from PCA, plot vector and color by tissue.
  # Maybe give fancy legend
  plot(y, main=title, col=rep(1:12, each=24), type='o') 
}

plot_periodogram <- function(x, N, title="Plot title"){
  # x is vector of data you want to check for rhythmicity
  N <- length(x)
  FF <- abs(fft(x) / sqrt(N)) ^ 2  # do my fft
  P <- FF[1:(N / 2 + 1)]  # only need first (N/2) + 1 values of FFT result
  f <- (0:(N / 2) / N)  # creates harmonic frequencies from 0 to 0.5 in steps of 1/12
  plot(f, P, type='l', main=title)
}



# MAIN --------------------------------------------------------------------


# define dirs
data_dir <- "microarray_data"
# fname <- "GSE54650_series_matrix.no_headers.txt"    # data straight from hogenesch
fname <- "hogenesch_2014_microarray.reprocessed.txt"    # data reprocessed by Affy package
fname <- "hogenesch_2014_rma.txt"    # data reprocessed by RMA package

# load data
data_path <- file.path(data_dir, fname)
dat <- read.table(data_path)

# log2 transform
dat_original <- dat
dat <- log2(dat)

# histogram expressions
dat_melt <- melt(dat)

ggplot(data=dat_melt, aes(x=value)) +
  geom_density()


# Plot quantiles ----------------------------------------------------------

p <- ggplot(dat_melt, aes(factor(variable), value))
p + geom_boxplot()


# Calculate PCA and Screeplot ---------------------------------------------


dat_standardize <- as.data.frame(scale(dat))

dat_pca <- prcomp(t(dat_standardize))

screeplot(dat_pca, type="lines", npcs = min(287, length(dat_pca$sdev)), log="y")


# Plot PCA ----------------------------------------------------------------


# Loop to plot several PCA components in series

for (x_comp in 1:5) {

  y_comp <- x_comp + 1
  # plot(dat_pca$x[, x_comp], dat_pca$x[, y_comp], main=paste("Component", x_comp, "vs.", y_comp))
  
  
  # Show tissue names only --------------------------------------------------
  
  
  newnames <- strsplit(names(dat_pca$x[,1]), '_')
  newnames <- sapply(newnames, function(x) {
    sample_id <- x[2]    # sample name with time point
    
    # THREE WAYS OF SHOWING SAMPLES NAMES:
    
    # 1.
    # show sample name with timepoint e.g. Liver34
    # tissue_name <- sample_id
    
    # 2.
    # Optional: remove time component of tissue name e.g. Liver
    tissue_name <- substring(sample_id, 1, nchar(sample_id) - 2)
    
    # 3.
    # Optional: show only time component e.g. 34
    # tissue_name <- substring(sample_id, nchar(sample_id) - 1, nchar(sample_id))
    })
  
  
  # Show text of tissues ----------------------------------------------------

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


# Loading plot ------------------------------------------------------------

# plot first n loadings
n <- 20
for (vector in 1:n){
  plot(dat_pca$x[, vector], main=paste("Vector", vector), col=rep(1:12, each=24), type='o')  
}



# Which vector loadings have oscillating components? ----------------------

N <- 288  # number of samples
T <- 24  # 24 hours in a period
# Make nicer rownames

rownames(dat_pca$x) <- newnames

# Fit PCA component of interest
for (pca_vector in 1:30){
  pca_component_string <- paste("PCA component:", pca_vector)
  # Create response vector, which is loadings
  
  y <- dat_pca$x[, pca_vector]
  
  # Create time vector, ignore different tissues
  # from t=0 to t=46 every 2 hours. Repeat for 12 tissues
  # t <- 1:288
  
  spec.pgram(y, log="no", main=pca_component_string)
  plot_loadings(y, title=paste("PCA component:", pca_vector))
  
  # t <- rep(seq(0, 46, 2), 12)
  t <- seq(0, 2 * N - 1, 2)
  omega <- (2 * pi) / (T)
  fit <- lm(y ~ cos(omega * t) + sin(omega * t))
  print("*********************************")
  print(paste0(pca_component_string))
  cat("\n")
  print(anova(fit))
  cat("\n")
}

pca_vector = 14
y <- dat_pca$x[, pca_vector]
FF <- abs(fft(y) / sqrt(N)) ^ 2
P <- FF[1:(N / 2 + 1)]
f <- (0:(N / 2) / N)
plot(f, P, type='l', main=paste("PCA Component:", pca_vector))







