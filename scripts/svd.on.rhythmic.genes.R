# Jake Yeung
# January 7 2015
# To eliminate genes with "noise", which shows up in projection and subsequent svd analysis, let's do projection only
# on genes that show rhythmicity in at least one tissue.
# 
# Dependencies: output of find.oscillating.genes.R is used to detect rhythmic genes.

library(plyr)
library(reshape2)
library(doMC)
library(PhaseHSV)
library(ggplot2)


# Parameters to play ------------------------------------------------------

## Specify the number of cores
mc <- 30

# Minimum log2 expression to be considered "expressed" in RNA-Seq
cutoff <- 5


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
source(file.path(scripts.dir, funcs.dir, "GetTopNValues.R"))
source(file.path(scripts.dir, funcs.dir, "OuterComplex.R"))
source(file.path(scripts.dir, funcs.dir, "OrderPhaseMatrix.R"))
source(file.path(scripts.dir, funcs.dir, "MergeToLong.R"))

Transform <- function(df, my.omega, normalize = TRUE){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # df: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
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
  
  transforms <- sapply(omegas, DoFourier, exprs = df$exprs, time = df$time)
  my.transformed <- transforms[which(omegas == omega)]  # corresponds to omega of interest
  
  if (normalize){
    # Normalize across omega
    # if median is 0, then set factor to 0, otherwise 1
    factor <- 1
    cutoff <- 5
    jmedian <- median(subset(df, experiment == "rnaseq")$exprs)
    if (jmedian <= cutoff){
      factor <- 0
    }
    
    my.transformed <- (my.transformed / sqrt(sum(Mod(transforms) ^ 2))) * factor
  } 
  return(data.frame(exprs.transformed = my.transformed))
}

GetInterval <- function(time.vector){
  # Given vector of times (equally spaced time points), return the time interval
  time.sort <- sort(unique(df$time))
  interval <- time.sort[2] - time.sort[1]
  return(interval)
}

IsRhythmic <- function(df, my.omega, pval.cutoff = 5e-3, method = "ANOVA"){
  # Test if rhythmic by BIC model selection: fit through all omegas.
  # df: long format, gene and condition
  # my.omega: omega of interest
  # method = "BIC" or "ftest"
  # pval.cutoff: for method = "ftest"
  
  fit.rhyth <- lm(exprs ~ 0 + experiment + sin(my.omega * time) + cos(my.omega * time), data = df)
  fit.flat <- lm(exprs ~ 0 + experiment, data = df)  # intercept only
  
  if (method == "BIC"){
    omegas.all <- GetOmegas(remove.zero = TRUE)
    # remove also my.omega to get omegas representing "noise"
    omegas.noise <- omegas.all[which(omegas.all != my.omega)]
    
    bic.test <- BIC(fit.rhyth, fit.flat)
    
    chosen.model <- rownames(bic.test)[which(bic.test$BIC == min(bic.test$BIC))]
    
    if (chosen.model == "fit.flat"){
      # if flat, no need to check for noise.
      return(FALSE)
    }
    # chosen model is fit.rhyth, but is it a noisy gene?
    # Check if noise components have an even better fit than fit.rhyth
    rhyth.bic <- bic.test["fit.rhyth", "BIC"]
    # fit noise
    fits.noise <- lapply(omegas.noise, 
                         function(w) lm(exprs ~ sin(w * time) + cos(w * time), 
                                        data = df))
    bic.noise.min <- min(sapply(fits.noise, BIC))
    # print(paste(unique(df$gene), "noise bic:", bic.noise.min, "rhyth bic:", rhyth.bic))
    if (rhyth.bic < bic.noise.min){
      return(TRUE)
    } else {
      return(FALSE)
    } 
  } else if (method == "ANOVA"){
    f.test <- anova(fit.flat, fit.rhyth)
    pval <- f.test[["Pr(>F)"]][[2]]
    if (is.nan(pval)){
      # if gene is all flat, then pval is NaN, force it to be 1 in this case.
      pval <- 1
    }
    if (pval < pval.cutoff){
      return(TRUE)
    } else {
      return(FALSE)
    }
}
}

ProjectToFrequency <- function(df, my.omega, normalize = TRUE, rhythmic.only = FALSE, method = "ANOVA", pval.cutoff = 5e-3){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # df: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
  # rhytmic.only: transforms only genes that are rhythmic (by BIC method)
  # 
  # omega in which we are interested.
  
  if (rhythmic.only){
    if (IsRhythmic(df, my.omega, pval.cutoff = pval.cutoff, method = "ANOVA")){
      df.transformed <- Transform(df, my.omega, normalize)
    } else {
      df.transformed <- data.frame(exprs.transformed = 0)
    }
  } else {
    df.transformed <- Transform(df, my.omega, normalize)
  }
  return(df.transformed)
}

# fit.rhyth <- lm(exprs ~ sin(omega * time) + cos(omega * time), data = df)
# fit.flat <- lm(exprs ~ 1, data = df)  # intercept only
# # omega.noise <- 2 * pi / (GetInterval(df$time) * 4)  # period of every two intervals is noise.
# omega.noise <- 2 * pi / 12
# fit.noise <- lm(exprs ~ sin(omega.noise * time) + cos(omega.noise * time), data = df)
# bic.test <- BIC(fit.rhyth, fit.flat, fit.noise)
# (bic.test)
# # extract rowname of minimum BIC
# chosen.model <- rownames(bic.test)[which(bic.test$BIC == min(bic.test$BIC))]
# if (chosen.model == "fit.rhyth"){
#   df.transformed <- Transform(df, my.omega, normalize)
# } else {
#   # if not rhythmic, then do not do transform it.
#   df.transformed <- data.frame(exprs.transformed = 0)
# }
# } else {
#   df.transformed <- Transform(df, my.omega, normalize)
# }

DoFourier <- function(exprs, time, omega, normalize = TRUE){
  p <- exprs %*% exp(-1i * omega * time)
  if (normalize){
    p <- p / length(time)  # normalize by number of datapoints 
  }
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


# Remove genes from RNA-Seqs that are not expressing. ---------------------

# # plot and show cutoff visually
# plot(density(unlist(rna.seq.exprs)))
# abline(v = cutoff, col = 'red')
# 
# jmean <- apply(rna.seq.exprs, 1, mean)
# 
# genes.not.expressed <- names(jmean[which(jmean < cutoff)])
# all.genes <- rownames(rna.seq.exprs)
# genes.expressed <- all.genes[which(!all.genes %in% genes.not.expressed)]
# 
# print(paste0(length(genes.not.expressed), 
#              " not expressed at cutoff = mean(", cutoff, 
#              ") for a gene across all tissues. Removing..."))
# 
# # Remove unexpressed genes
# 
# rna.seq.exprs.filtered <- rna.seq.exprs[genes.expressed, ]
rna.seq.exprs.filtered <- rna.seq.exprs  # no filter

# Take only common genes --------------------------------------------------

array.genes <- rownames(normalized.array)
rna.seq.genes <- rownames(rna.seq.exprs.filtered)

common.genes <- intersect(array.genes, rna.seq.genes)

print(paste(length(common.genes), "common genes between array and rnaseq"))

normalized.array <- normalized.array[common.genes, ]
rna.seq.exprs.filtered <- rna.seq.exprs.filtered[common.genes, ]

# Merge data into long format ---------------------------------------------

dat <- MergeToLong(normalized.array, rna.seq.exprs.filtered)

dat.split <- split(dat, dat$tissue)
 

# Project my data ---------------------------------------------------------

omega <- 2 * pi / 24
omegas <- GetOmegas()

# # # test my function
# problem genes
problem.genes <- c("Hephl1", "Dnmt3l", "Gm10804", "Fgf14", "Spint4", "Defb23", "Defb34", "B3galt2")
problem.genes <- c("B3galt2", "Sgcg", "Rgs16", "Ddit4", "Smpx", "Eef1a2", "D3Bwg0562e", "Slc9a2", 
                   "Tceal3", "Murc", "Fgf10", "Dclk1", "Ccne1", "Cox6a2", "Trim63", 
                   "Kcng4", "Slc17a9", "Diras2", "Txlnb", "A830018L16Rik", "Elovl3", "Alb")

problem.genes <- c("Diras2", "Kcng4", "Trim63", "Cox6a2", "Alb", "Elovl3", "Spint4", "Svs5")

df <- subset(dat.split$WFAT, gene %in% c(clockgenes, problem.genes))
df.proj <- ddply(df, .(gene), ProjectToFrequency, omega, normalize = FALSE, rhythmic.only = TRUE)
df.proj$mod <- Mod(df.proj$exprs.transformed)
df.proj <- df.proj[order(df.proj$mod, decreasing = TRUE), ]
(df.proj)
# 
# # test individual genes
# jgene <- "Gm10233"
# jgene <- "Mir181b.2"
# jgene <- "Gje1"
jgene <- "B3galt2"
jgene <- "Arntl"
jmean <- mean(subset(dat.split$Liver, gene == jgene & experiment == "rnaseq")$exprs)
transforms <- sapply(omegas, DoFourier, exprs = subset(dat.split$Liver, gene == jgene)$exprs, time = subset(dat.split$Liver, gene == jgene)$time)


start.time <- Sys.time()
if (getDoParWorkers() == 1){
  registerDoMC(40)
}
print(paste("Parallel processing on:", getDoParWorkers(), "cores"))
dat.split.proj <- lapply(dat.split, function(x){
  ddply(x, .(gene), ProjectToFrequency, my.omega = omega, normalize = FALSE, rhythmic.only = TRUE, pval.cutoff = 5e-5, .parallel = TRUE)
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
rownames(s$v) <- colnames(dat.wide)

# Plot interesting components ---------------------------------------------

tissues <- GetTissues(colnames(normalized.array))

sing.vals <- seq(length(s$d))
# sing.vals <- c(1, 2, 3)

pdf("plots/components.rhythmic.only.filter.noise.ANOVA.5e5.rotate.pdf")
# ScreePlot
plot(s$d^2 / sum(s$d ^ 2), type='o')  # Manual screeplot. 

for (sing.val in sing.vals){
  print(paste("Plotting component:", sing.val))
  eigengene <- s$v[, sing.val]
  jmax <- max(Mod(eigengene))
  max.loading <- names(eigengene[which(Mod(eigengene) == jmax)])
  jmax.arg <- Arg(eigengene[max.loading])  # aligns everything relative to jmax.arg
  PlotComplex(eigengene, 
              axis.min = -jmax,
              axis.max = jmax,
              labels = tissues, 
              col = "HSV",
              add.text.plot = FALSE, 
              main = paste("Component:", sing.val),
              jpch = 20,
              threshold = 0, 
              rotate = -jmax.arg)
  
  eigensample <- s$u[, sing.val]
  jmax <- max(Mod(eigensample))
  # rotate opposite way to make it all kosher
  PlotComplex(eigensample,
              axis.min = -jmax,
              axis.max = jmax,
              labels = names(eigensample),
              col = "HSV",
              main = paste("Component:", sing.val),
              add.text.plot = FALSE,
              jpch = 1,
              threshold = 0.5 * jmax,
              verbose = FALSE,
              rotate = jmax.arg)
  
  # show how fast the values drop
  top.genes.all <- GetTopNValues(Mod(eigensample), N = length(eigensample))
  plot(top.genes.all$vals, 
       log = "x", 
       main = paste("Component:", sing.val),
       ylab = "Modulus of gene loading",
       xlab = "Gene index")
  # Plot args in color
  top.genes <- GetTopNValues(Mod(eigensample), N = 100)# list of $vals $i
  # use drop to keep rownames
  outer.prod.mat <- s$d[sing.val] * OuterComplex(s$u[top.genes$i, sing.val, drop = FALSE], t(s$v[, sing.val, drop = FALSE]))
  outer.prod.mat <- OrderPhaseMatrix(outer.prod.mat, order.by = max.loading, order.tissues = TRUE)
  PlotArgsMatrix(outer.prod.mat, main = paste("Component", sing.val))
  
  # print top 100 genes (copy and paste-able)
  cat(paste0(head(names(top.genes$vals), n = 50), "\n"))
  
#   # print top 20 genes (list-able)
  cat(paste(head(names(top.genes$vals), n = 20), collapse = '", "'))
  cat("\n")
}
dev.off()

# plot rhythmic genes
pdf("plots/rhythmic.genes.by.singular.values.rhythmic.only.filter.noise.ANOVA.5e5.pdf")
for (sing.val in sing.vals){
  eigensample <- s$u[, sing.val]
  top.genes <- GetTopNValues(Mod(eigensample), N = 10)# list of $vals $i
  
  genes <- names(top.genes$vals)
  
  dat.sub <- subset(dat, gene %in% genes)
  dat.sub$gene <- factor(dat.sub$gene)
  
  for (jgene in names(top.genes$vals)){
    print(jgene)
    print(PlotGeneAcrossTissues(subset(dat.sub, gene == jgene), jtitle = paste("gene:", jgene, "component", sing.val)))
  }
}
dev.off()