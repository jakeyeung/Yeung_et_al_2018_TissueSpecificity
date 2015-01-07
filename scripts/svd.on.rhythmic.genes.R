# Jake Yeung
# January 7 2015
# To eliminate genes with "noise", which shows up in projection and subsequent svd analysis, let's do projection only
# on genes that show rhythmicity in at least one tissue.
# 
# Dependencies: output of find.oscillating.genes.R is used to detect rhythmic genes.

library(ddply)
library(reshape2)
library(doMC)

## Specify the number of cores
registerDoMC(30)

## Check how many cores we are using
getDoParWorkers()

clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')

# My Functions ------------------------------------------------------------

scripts.dir <- "scripts"
funcs.dir <- "functions"
source(file.path(scripts.dir, funcs.dir, "MergeAndProject.R"))
source(file.path(scripts.dir, funcs.dir, "GetTissueTimes.R"))
source(file.path(scripts.dir, funcs.dir, "FourierFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "LoadAndHandleData.R"))
source(file.path(scripts.dir, funcs.dir, "DataHandlingFunctions.R"))
source(file.path(scripts.dir, funcs.dir, "PlotFunctions.R"))

ProjectToFrequency <- function(df, my.omega, normalize = TRUE){
  # Perform fourier transform and normalize across all frequencies.
  # 
  # Input:
  # df: long dataframe containing expression and time columns
  # to be transformed to frequency domain.
  # 
  # omega in which we are interested.
  
  # if we normalize, get list of omegas.
  # otherwise we just use omega
  if (normalize){
    # t <- sort(unique(df$time))
    # n.timepoints <- length(t)
    # interval <- t[2] - t[1]
    omegas <- GetOmegas(remove.zero = TRUE)
  } else {
    omegas <- my.omega
  }
  
  sum.transforms <- 0
  
  for (omega in omegas){
    # take gene from a tissue and get its fourier transform
    transformed <- df$exprs %*% exp(-1i * omega * df$time)
    # normalize by number of datapoints
    transformed <- transformed / length(df$exprs)
    sum.transforms <- sum.transforms + Mod(transformed)
    
    print(transformed)
    
    if (omega == my.omega){
      my.transformed <- transformed
    }
  }
  
  if (normalize){
    # Normalize across omegas
    my.transformed <- my.transformed / sum.transforms 
  } 
  
  data.frame(exprs.transformed = my.transformed)
}

DoFourier <- function(exprs, time, omega){
  p <- exprs %*% exp(-1i * omega * time)
  return(p)
}

GetOmegas <- function(n.timepoints=24, interval=2, remove.zero=FALSE){
  # Get omegas from a time series
  # 
  # remove.zero: removes first element (omega = 0).
  # Useful for normalizing Fourier when I don't want omega = 0 case
  
  # Output:
  # list of omegas: can be used for fourier transform
  
  t <- seq(from = 0, by = 1, length.out = n.timepoints)  # will account for interval later
  
  # get harmonic frequencies
  # this creates harmonic frequencies from 0 to .5 in steps of 1/n.timepoints
  
  mid.point <- length(t) / 2 + 1
  freqs <- t[1:mid.point] / n.timepoints
  
  # Convert to period before adjusting for interval
  T.unadj <- 1 / freqs
  
  # Account for interval
  T <- T.unadj * interval
  
  # Calculate omegas
  omegas <- 2 * pi / T
  
  if (remove.zero){
    omegas <- omegas[2:length(omegas)]
  }
  
  return(omegas)
}

ConvertLongToWide <- function(long.df, measurement.var = "exprs.transformed"){
  wide.df <- dcast(dat.proj, gene ~ tissue, value.var = measurement.var)
  # first row is gene name, let's make them rowname and remove first column.
  
  rownames(wide.df) <- wide.df$gene
  
  wide.df <- subset(wide.df, select = -c(gene))
  
  return(wide.df)
}

PlotGeneAcrossTissues <- function(dat, jtitle){
  m <- ggplot(dat, aes(x = time, y = exprs,
                   group = experiment, 
                   colour = experiment)) + 
      geom_point() + 
      geom_line() + 
      facet_wrap(~tissue) +
      ggtitle(jtitle)
  return(m)
}



# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
normalized.array.fname <- "array.adj.0.07.txt"
normalized.array.path <- file.path(data.dir, normalized.array.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)

# Load file ---------------------------------------------------------------

normalized.array <- LoadNormalizedArray(normalized.array.path)
rna.seq.exprs <- LoadRnaSeq(rna.seq.path)


# Log2 transform of array and rnaseq --------------------------------------

normalized.array <- log2(normalized.array + 1)
rna.seq.exprs <- log2(rna.seq.exprs + 1)


# Take only common genes --------------------------------------------------

array.genes <- rownames(normalized.array)
rna.seq.genes <- rownames(rna.seq.exprs)

common.genes <- intersect(array.genes, rna.seq.genes)

print(paste(length(common.genes), "common genes between array and rnaseq"))

normalized.array <- normalized.array[common.genes, ]
rna.seq.exprs <- rna.seq.exprs[common.genes, ]

# Merge data into long format ---------------------------------------------

dat <- MergeToLong(normalized.array, rna.seq.exprs)

dat.split <- split(dat, dat$tissue)


# Project my data ---------------------------------------------------------

omega <- 2 * pi / 24

# # test my function
clockgenes <- c(clockgenes, "Spint4", "Defb23")
df <- subset(dat.split$WFAT, gene %in% clockgenes)
df.proj <- ddply(df, .(gene), ProjectToFrequency, omega, normalize = TRUE)

start.time <- Sys.time()
dat.split.proj <- lapply(dat.split, function(x){
  ddply(x, .(gene), ProjectToFrequency, my.omega = omega, normalize = TRUE, .parallel = TRUE)
})
print(Sys.time() - start.time)


# Add tissue information into each list -----------------------------------

for (tissue in names(dat.split.proj)){
  dat.split.proj[[tissue]]$tissue <- tissue
}

# Combine data ------------------------------------------------------------

dat.proj <- do.call(rbind, dat.split.proj)


# Remove NaNs -------------------------------------------------------------

# dat.proj <- dat.proj[is.finite(dat.proj$exprs.transformed), ]

# Clean up ----------------------------------------------------------------

# rm(dat.split.proj, dat.split)

# Do SVD on this normalized dataset ---------------------------------------

# long to wide conversion
dat.wide <- ConvertLongToWide(dat.proj, measurement.var = "exprs.transformed")

# Complete cases only
dat.wide <- dat.wide[complete.cases(dat.wide), ]

s <- svd(dat.wide)


# Set rownames for u and v ------------------------------------------------

rownames(s$u) <- rownames(dat.wide)
colnames(s$v) <- colnames(dat.wide)

# Plot interesting components ---------------------------------------------

sing.vals <- seq(length(s$d))

# ScreePlot
plot(s$d^2, type='o')  # Manual screeplot. 

for (sing.val in sing.vals){
  eigengene <- s$v[, sing.val]
  jmax <- max(Mod(eigengene))
  PlotComplex(eigengene, 
              axis.min = -jmax,
              axis.max = jmax,
              labels = tissues, 
              col = "HSV",
              add.text.plot = FALSE, 
              main = paste("Component:", sing.val),
              jpch = 20,
              threshold = 0)
  
  eigensample <- s$u[, sing.val]
  jmax <- max(Mod(eigensample))
  PlotComplex(eigensample,
              axis.min = -jmax,
              axis.max = jmax,
              labels = names(eigensample),
              col = "HSV",
              main = paste("Component:", sing.val),
              add.text.plot = FALSE,
              jpch = 1,
              threshold = 0.75 * jmax,
              verbose = TRUE)
}

# genes <- c("Dbp", "Defb22", "Defb23", "Defb34", "Npas2", "Nr1d1", "Olfr1014", "Olfr1015", "Olfr1052", "Olfr731", "Spint4")
# 
# dat.sub <- subset(dat, gene %in% genes)
# for (jgene in genes){
#   m <- PlotGeneAcrossTissues(subset(dat.sub, gene == jgene), jtitle = jgene)
#   print(m)
# }



