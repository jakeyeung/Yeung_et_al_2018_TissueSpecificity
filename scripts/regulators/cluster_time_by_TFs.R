# cluster_time_by_TFs.R
# Jake Yeung
# March 12 2015

library(ggplot2)
library(plyr)
library(reshape2)
library(gplots)
library(doMC)

# Paths -------------------------------------------------------------------

# define dirs
data.dir <- "data"
tf.fname <- "motifs_and_TFs.list"
tf.path <- file.path(data.dir, tf.fname)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/SvdFunctions.R")
source('scripts/functions/PlotGeneAcrossTissues.R')

# Main --------------------------------------------------------------------

dat <- LoadArrayRnaSeq()
tfs <- GetTFs()

dat.tfs <- subset(dat, gene %in% tfs)

dat.split <- split(dat.tfs, dat.tfs$tissue)

dat.split$WFAT <- NULL

omega <- 2 * pi / 24
omegas <- GetOmegas()

start.time <- Sys.time()
if (getDoParWorkers() == 1){
  registerDoMC(40)
}
print(paste("Parallel processing on:", getDoParWorkers(), "cores"))
dat.split.proj <- lapply(dat.split, function(x){
  ddply(x, .(gene), ProjectToFrequency, my.omega = omega, normalize = FALSE, rhythmic.only = FALSE, pval.cutoff = 1, .parallel = TRUE)
})
print(Sys.time() - start.time)

# Add tissue information into each list -----------------------------------

for (tissue in names(dat.split.proj)){
  dat.split.proj[[tissue]]$tissue <- tissue
}

# Combine data ------------------------------------------------------------

dat.proj <- do.call(rbind, dat.split.proj)

dat.proj$mod <- Mod(dat.proj$exprs.transformed)

# ggplot(subset(dat.proj, !(tissue %in% c("BS", "Cere", "Hypo"))), aes(x = mod)) + geom_density() + xlim(0, 0.1)
ggplot(subset(dat.proj, mod > 0), aes(x = mod)) + geom_density()

# long to wide conversion
dat.wide <- ConvertLongToWide(dat.proj, measurement.var = "mod")

# Complete cases only. This removes NaN rows.
dat.wide <- dat.wide[complete.cases(dat.wide), ]


# Heatmap -----------------------------------------------------------------
my.palette <- colorRampPalette(c("black", "yellow"))(n = 299)

# # (optional) defines the color breaks manually for a "skewed" color transition
col.breaks = c(seq(0, 5, length=150),  # black
               seq(5, 20, length=151))  # yellow

dat.wide.sub <- dat.wide[c("Arntl", "Nr4a2", "Ppara", "Hnf4a", "Myod1", "Snai3", "Id1", "Stat4"), ]
heatmap.2(as.matrix(dat.wide.sub),
          # cellnote = dat.tfs.mat,  # same data set for cell labels
          main = "Rhythmic exprs", # heat map title
          # notecol="black",      # change font color of cell labels to black
          density.info='density',  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5, 45),     # widens margins around plot
          dendrogram = "both",
          cexRow = 1,
          col=my.palette,       # use on color palette defined earlier 
          # breaks=col.breaks    # enable color transition at specified limits
          # dendrogram="col",     # only draw a row dendrogram
          # Colv="NA")            # turn off column clustering
)
