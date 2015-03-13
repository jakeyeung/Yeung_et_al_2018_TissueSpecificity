# cluster_tissues_by_TFs.R
# Jake Yeung
# March 9 2015

library(ggplot2)
library(plyr)
library(reshape2)
library(gplots)

# Paths -------------------------------------------------------------------

# define dirs
data.dir <- "data"
tf.fname <- "motifs_and_TFs.list"
tf.path <- file.path(data.dir, tf.fname)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/FirstColToRowname.R")

AvgExprsAcrossTissues <- function(df){
  # take average exprs from df which contains exprs colname 
  # df.avg.tiss <- ddply(df, .(tissue), AvgExprs)
  df.avg.tiss <- ddply(df, .(tissue), summarise,
                       mean.exprs = mean(exprs))
  return(df.avg.tiss)
}

AllOnes <- function(v){
  # check if vector is all 1s
  return(all(v == 1))
}

FilterAllOnesOrZeros <- function(mat){
  # Filter out a binary matrix only for rows that do not contain
  # all 1s or 0s. They are not meaningful for differentiating TF
  # across tissues
  mat <- mat[which(apply(mat, 1, AllOnes) == FALSE), ]  # remove TFs with only 1s 
  mat <- mat[which(rowSums(mat) > 0), ] # remove TFs with only 0s
  return(mat)
}

Binarize <- function(mat, cutoff){
  mat[which(mat < cutoff)] <- 0
  mat[which(mat >= cutoff)] <- 1
  return(mat)
}


# Set pallette for heatmap colors -----------------------------------------


# creates a own color palette from red to green
# my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
my.palette <- colorRampPalette(c("black", "yellow"))(n = 299)

# # (optional) defines the color breaks manually for a "skewed" color transition
col.breaks = c(seq(0, 5, length=150),  # black
               seq(5, 20, length=151))  # yellow


# Main --------------------------------------------------------------------

dat <- LoadArrayRnaSeq()

# plot distribution of expression 
ggplot(aes(x = exprs, fill = experiment), data = dat) +
  geom_density(alpha = 0.3)

tfs <- GetTFs()

dat.tfs <- subset(dat, gene %in% tfs & experiment == "rnaseq")

dat.tfs.avg <- ddply(dat.tfs, .(gene), AvgExprsAcrossTissues)

dat.tfs.mat <- dcast(dat.tfs.avg, gene ~ tissue, value.var = "mean.exprs")

dat.tfs.mat <- as.matrix(FirstColToRowname(dat.tfs.mat))

# binarize
cutoff <- 5
dat.tfs.mat.binary <- Binarize(dat.tfs.mat, cutoff)
dat.tfs.mat.binary <- FilterAllOnesOrZeros(dat.tfs.mat.binary)


# heatmap that shiet
heatmap.2(dat.tfs.mat.binary,
          # cellnote = dat.tfs.mat,  # same data set for cell labels
          main = paste0("exprs cutoff=", cutoff), # heat map title
          # notecol="black",      # change font color of cell labels to black
          density.info='density',  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          # margins =c(5, 45),     # widens margins around plot
          dendrogram = "both",
          cexRow = 0.3,
          col=my.palette,       # use on color palette defined earlier 
          # breaks=col.breaks    # enable color transition at specified limits
          # dendrogram="col",     # only draw a row dendrogram
          # Colv="NA")            # turn off column clustering
)
# heatmap that shiet
heatmap.2(dat.tfs.mat,
          # cellnote = dat.tfs.mat,  # same data set for cell labels
          main = "TF exprs", # heat map title
          # notecol="black",      # change font color of cell labels to black
          density.info='density',  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5, 45),     # widens margins around plot
          dendrogram = "both",
          cexRow = 0.1,
          col=my.palette,       # use on color palette defined earlier 
          # breaks=col.breaks    # enable color transition at specified limits
          # dendrogram="col",     # only draw a row dendrogram
          # Colv="NA")            # turn off column clustering
)