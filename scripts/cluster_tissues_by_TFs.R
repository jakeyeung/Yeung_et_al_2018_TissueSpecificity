# cluster_tissues_by_TFs.R
# Jake Yeung
# March 9 2015

library(ggplot2)
library(plyr)
# Paths -------------------------------------------------------------------

# define dirs
data.dir <- "data"
tf.fname <- "motifs_and_TFs.list"
tf.path <- file.path(data.dir, tf.fname)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")

AvgExprs <- function(df){
  return(mean(df$exprs))  
}

AvgExprsAcrossTissues <- function(df){
  # take average exprs from df which contains exprs colname 
  # df.avg.tiss <- ddply(df, .(tissue), AvgExprs)
  df.avg.tiss <- ddply(df, .(tissue), summarise,
                       mean = mean(exprs))
  return(df.avg.tiss)
}

# Main --------------------------------------------------------------------

dat <- LoadArrayRnaSeq()

# plot distribution of expression 
ggplot(aes(x = exprs, fill = experiment), data = dat) +
  geom_density(alpha = 0.3)

tfs <- GetTFs()

dat.tfs <- subset(dat, gene %in% tfs & experiment == "rnaseq")

dat.tfs.avg <- ddply(dat.tfs, .(gene), AvgExprsAcrossTissues)

