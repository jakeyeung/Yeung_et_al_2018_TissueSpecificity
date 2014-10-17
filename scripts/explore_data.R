# Jake Yeung
# explore_data.R
# initial exploration of data in microarray across tissues

library(ggplot2)
library(reshape)
library(wordcloud)  # for showing text without jumbling

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

for (vector in 1:20){
  plot(dat_pca$x[, vector], main=paste("Vector", vector), col=rep(1:12, each=24))  
}


