# Jake Yeung
# explore_data.R
# initial exploration of data in microarray across tissues

library(ggplot2)
library(reshape)
library(wordcloud)  # for showing text without jumbling


# Define my functions -----------------------------------------------------

# First source my functions I wrote
funcs.dir <- file.path('scripts', 'functions')
source(file.path(funcs.dir, 'SampleNamehandler.R'))  # for shortening sample names
source(file.path(funcs.dir, 'PcaPlotFunctions.R'))  # for visualizing PCA, periodoigrams
source(file.path(funcs.dir, 'FourierFunctions.R'))  # for periodoigrams

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
data_dir <- "data"
# fname <- "hogenesch_2014_rma.txt"    # data reprocessed by RMA package
fname <- "hogenesch_2014_rma.genenames.colnameordered.txt"  # new with gene names

# load data
data_path <- file.path(data_dir, fname)
print(paste("Reading data from,", data_path, "May take a few a minutes."))
dat <- read.table(data_path, header=TRUE, sep='\t')
print("Read data to memory.")

# Handle duplicate rownames: MICROARRAY -------------------------------------

# first column contained gene names
rownames(dat) <- make.names(dat$gene, unique=TRUE)
# genes <- dat$gene
# coords <- dat$coordinates
# probeid <- dat$ID_REF

drop.cols <- c("gene")
dat <- dat[, !(names(dat) %in% drop.cols)]
# Peek(dat)


# Get colnames ------------------------------------------------------------

# make more meaningful sample names
# user changeable paramters: show="tissue.time" | "tissue" | "time"
colnames(dat) <- ShortenSampNames(colnames(dat), show="tissue")

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

dat_pca <- prcomp(t(dat), center=TRUE, scale.=TRUE)

screeplot(dat_pca, type="lines", npcs = min(287, length(dat_pca$sdev)), log="y")


# Garbage collection ------------------------------------------------------

rm(dat)


# Plot PCA: tissue components ----------------------------------------------------

# Loop to plot scatter plot of PCA i versus PCA i + 1
for (x_comp in 1:3) {

  y_comp <- x_comp + 1  # plot x PCA component i on x-axis, PCA component i+1 on y-axis
  

  # Show PCA plot of tissues ----------------------------------------------------
  
  # TWO WAYS OF COLORING SAMPLES IN PCA: by time or by tissue
  
  
  # 1.
  #   # Color by time point
  #   colors.by.time <- rep(1:24, 12)
  #   textplot(dat_pca$x[, x_comp], dat_pca$x[, y_comp], rownames(dat_pca$x), cex=0.7, col=colors.by.time, 
  #           main=paste("Component", x_comp, "vs.", y_comp))
  
  # 2. 
  # Color by tissue
  colors.by.tissue <- rep(1:12, each=24)
  textplot(dat_pca$x[, x_comp], dat_pca$x[, y_comp], rownames(dat_pca$x), cex=0.7, col=colors.by.tissue, 
           main=paste("Component", x_comp, "vs.", y_comp))
  
}


# Plot PCA: circadian components ---------------------------------------------


# Loop to plot scatter plot of PCA i versus PCA i + 1
# FELIX: this is for the 'circadian; components
for (x_comp in 17:19) {
  
  y_comp <- x_comp + 1  # plot x PCA component i on x-axis, PCA component i+1 on y-axis
  
  
  # Show PCA plot of tissues ----------------------------------------------------
  
  # TWO WAYS OF COLORING SAMPLES IN PCA
  
  # 1.
  #   # Color by time point
  colors=hsv(rep(c(1:12)/12, 24), 1, 1)

  ii=1:nrow(dat_pca$x)
#   ii=which(newnames=="Liver")

  newnames.2=paste(colnames(dat_pca$x), rep(c(0:23)*2, 12))
  
  r=max(range(c(dat_pca$x[ii, x_comp], dat_pca$x[ii, y_comp]))) * 1.2

  textplot(dat_pca$x[ii, x_comp], dat_pca$x[ii, y_comp], newnames.2[ii], cex=0.7, col=colors, 
           main=paste("Component", x_comp, "vs.", y_comp), xlim=c(-r,r), ylim=c(-r,r))
}

# Which vector loadings have oscillating components? ----------------------
# We will fit two linear models:
# 1) linear model with common intercept and common oscillating component
# 2) linear model with tissue-specific intercept 
# and common oscillating component

# Optional: can MeanCenter the vector y to remove cross-tissue differences
# and focus only on oscillations.

N <- nrow(dat_pca$x)  # number of samples.
T <- 24  # 24 hours in a period

# Fit PCA component of interest: loop to try many different PCAs
# user changeable range
for (pca_vector in 10:20) {
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
  
  tissue.factors <- GetTissueSpecificMatrix(n.tissues, n.timepoints)
  
  # fit my lm using cos and sin with angular frequency
  fit <- lm(y ~ cos(omega * t) + sin(omega * t))
  fit.multi <- lm(y ~ 0 + tissue.factors + cos(omega * t) + sin(omega * t))
  # END: Linear fit for period of 24 hours
  
  # Print some statementsy to describe fit
  cat("*********************************\n")
  cat(paste0("PCA ", pca_vector, "\n"))
  cat("\n")
  
  cat("Max freqs\n")
  print(max.freqs)
  
  cat("Simple fit:\n")
  print(anova(fit))
  cat("\n")
  
  cat("Multi fit:\n")
  print(anova(fit.multi))
  cat("\n")
  
  cat("Comparing simple and multi fit:\n")
  print(anova(fit, fit.multi))
  cat("\n")
  
}


# Get tissue-specific genes -----------------------------------------------

PCA <- 2
gene.loadings <- dat_pca$rotation[, PCA]
# take top absolute values
top.genes <- order(gene.loadings)[1:10]
# print top genes
rownames(dat_pca$rotation)[top.genes]
