# 2015-06-17
# Jake Yeung
# adjust_array_to_kallisto.R
# Adjust array to kallisto gene counts because I trust it more than DESeq

setwd("/home/yeung/projects/tissue-specificity")
library(ggplot2)
library(dplyr)

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/RegressionFunctions.R")

saturation <- function(params, r) params[2] + (params[1] * r) / (params[3] + r)
saturation.inv <- function(params, m) params[3] * (m - params[2]) / (params[2] + params[1] - m)

# LINEAR FUNCTIONS
linear <- function(params, x) params[1] + params[2] * x
linear.inv <- function(params, y) {
  if (is.na(params[2])){
    params[2] <- Inf
  }
  (y - params[1]) / params[2]
}


GetGeneList <- function(){
  # for testing
  # Clock genes --------------------------------------------------------
  
  clockgenes <- c('Nr1d1','Dbp', 'Arntl', 'Npas2', 'Nr1d2', 
                  'Bhlhe41', 'Nfil3', 'Cdkn1a', 'Lonrf3', 
                  'Tef', 'Usp2', 'Wee1', 'Dtx4', 'Asb12')
  clockgenes <- c(clockgenes, 'Elovl3', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry2', 'Cry1')
  clockgenes <- c(clockgenes, 'Defb48', 'Svs1', 'Svs2', 'Svs5', 'Defb20', 'Adam7', 'Lcn8', 'Rnase10', 'Teddm1')
  clockgenes <- c(clockgenes, "Sry")  # known to not be expressed via RNASeq
  
  
  # Tissue genes ------------------------------------------------------------
  
  clockgenes <- c("Plbd1","Acacb","Hadh","Decr1","Acadl","Ech1","Acsl1","Cpt2","Acaa2")  # PCA 1 from explore_data.R
  clockgenes <- c(clockgenes, "Elovl1","Slc35b3","Fkbp15","Slc39a8","Sep15","Mospd2","Med11","Slco2a1","Arf6","Yipf6")  # PCA 2 from explore_data.R
  
  
  # Problem genes -----------------------------------------------------------
  
  # interesting cases
  clockgenes <- c(clockgenes, "Tmem59l", "Saysd1", "H2.Q1", "Col25a1", "Clec18a", "4930427A07Rik", 
                        "Ncan", "Crtam", "Fam43b", "Nphp4", "Nuf2", "Sox11", "Krt23", "Myo1h", 
                        "Mas1", "Cd207", "Tgif2", "Sdsl", "Gm8659", "Fsd1", "X2510049J12Rik", 
                        "1600029I14Rik", "Syndig1l", "Cyp4x1", "E130309D14Rik", "Gm281", 
                        "Wdr27", "Daw1", "Tnfsf9", "Myo16")
  
  # genes with 0 expression
  clockgenes <- c(clockgenes, "1700025F22Rik")
  
  return(clockgenes)
}

LoadKallistoGene <- function(inpath, gene_colname = "gene_name", log2.pseudocount=FALSE){
  source("scripts/functions/ConvertRNASeqTissueNamesToArray.R")
  
  
  if (missing(inpath)){
    inpath <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
  }
  
  dat <- read.table(inpath, header = TRUE)
  
  genes <- dat[, gene_colname]
  dat.mat <- dat[, colnames(dat)[which(colnames(dat) != gene_colname)]]
  
  if (log2.pseudocount != FALSE){
    dat.mat <- log2(dat.mat + log2.pseudocount)  # add a pseudocount
  }
  
  tissues <- sapply(colnames(dat.mat), function(s) strsplit(s, '_')[[1]][[1]])
  tissues <- ConvertRNASeqTissueNamesToArray(tissues)
  times <- GetTimes(colnames(dat.mat), get_unique=FALSE)
  
  tpm.long <- data.frame(gene = rep(genes, ncol(dat.mat)),
                         tissue = rep(tissues, each = nrow(dat.mat)),
                         time = as.numeric(rep(times, each = nrow(dat.mat))),
                         tpm = unlist(dat.mat))
}

PlotGeneTpm <- function(dat, jgene, log2.pseudocount=FALSE, scale=FALSE){
  dat.sub <- subset(dat, gene == jgene)
  if (scale != FALSE){
    dat.sub$tpm <- dat.sub$tpm * scale
  }
  if (log2.pseudocount != FALSE){
    dat.sub$tpm <- log2(dat.sub$tpm + log2.pseudocount)
  }

  ggplot(dat.sub, aes(x = time, y = tpm)) + geom_point() + geom_line() + facet_wrap(~tissue)
}

LoadArray <- function(inpath, gene_colname = "gene", get.norm = FALSE, form = "long"){
  if (missing(inpath)){
    inpath <- "/home/yeung/projects/tissue-specificity/data/array_exprs_colnames_fixed.best.probe.selected.txt"
  }
  dat <- read.table(inpath, header=TRUE)
  dat.mat <- dat[, colnames(dat)[which(colnames(dat) != gene_colname)]]
  if (form == "wide"){
    rownames(dat) <- dat[[gene_colname]]
    dat[[gene_colname]] <- NULL
  } 
  if (form == "wide" & get.norm == FALSE){
    return(dat)
  } else if(form == "wide" & get.norm == TRUE){
    # genes <- as.character(dat[[gene_colname]])
    # dat.mat <- 2 ^ dat.mat
    return(2 ^ dat)
  }  
  tissues <- GetTissues(colnames(dat.mat), get_unique = FALSE)
  times <- GetTimes(colnames(dat.mat), get_unique = FALSE)
  
  genes <- dat[[gene_colname]]
  array.long <- data.frame(gene = rep(genes, ncol(dat.mat)),
                           tissue = rep(tissues, each = nrow(dat.mat)),
                           time = as.numeric(rep(times, each = nrow(dat.mat))),
                           signal = unlist(dat.mat))
  if (get.norm == TRUE){
    array.long$signal.norm <- 2 ^ array.long$signal
  }
  return(array.long)
}


GetMeanVarByTissues <- function(exprs, tissue.names){
  # Calculate mean and variance for each gene per tissue
  N <- nrow(exprs) * length(tissue.names)  # one measurement for each gene for all tissues.
  mean.var <- list(mean=matrix(NA, nrow=nrow(exprs), ncol=length(tissue.names),
                               dimnames = list(rownames(exprs), 
                                               tissue.names)), 
                   var=matrix(NA, nrow=nrow(exprs), ncol=length(tissue.names),
                              dimnames = list(rownames(exprs), tissue.names)))
  for (j in 1:length(tissue.names)){
    tissue <- tissue.names[j]
    gene.tissue.exprs <- exprs[, grepl(tissue, colnames(exprs))]
    # calculate mean and var, by row
    exprs.mean <- apply(gene.tissue.exprs, 1, mean)
    exprs.var <- apply(gene.tissue.exprs, 1, var)
    # append to matrix 
    mean.var$mean[, j] <- exprs.mean
    mean.var$var[, j] <- exprs.var
  }
  # Make dataframe for ggplot2 and other analyses
  mean.var.df <- data.frame(mean=as.vector(mean.var$mean), var=as.vector(mean.var$var))
  return(mean.var.df)
}

PlotArrayRnaseq <- function(dat, log2.transform=FALSE){
  if (log2.transform  != FALSE){
    dat$signal.norm <- log2(dat$signal.norm + 1)
  }
  ggplot(dat, aes(x = time, y = signal.norm, colour = experiment, group = experiment, fill = experiment)) + 
    geom_point() +
    geom_line() +
    facet_wrap(~tissue)
}



MergeRnaseqArray <- function(kallisto.long, array.long, array.var = "signal.norm"){
#   kallisto.long$experiment <- rep("rnaseq", nrow(kallisto.long))
#   array.long$experiment <- rep("array", nrow(array.long))
  
  # Get common genes
  rnaseq.genes <- unique(kallisto.long$gene)
  array.genes <- unique(array.long$gene)
  common.genes <- intersect(rnaseq.genes, array.genes)
  
  # We only want the overlap of times, RNA-Seq is subset of Array times
  rnaseq.times <- unique(kallisto.long$time)
  array.times <- unique(array.long$time)
  common.times <- intersect(rnaseq.times, array.times)
  
  kallisto.sub <- subset(kallisto.long, gene %in% common.genes & time %in% common.times)
  array.sub <- subset(array.long, gene %in% common.genes & time %in% common.times)
  
  # Fix factor levels
  kallisto.sub$gene <- factor(kallisto.sub$gene)
  array.sub$gene <- factor(array.sub$gene)
  
  # order them up
  kallisto.sub <- kallisto.sub[with(kallisto.sub, order(gene, tissue, time)), ]
  array.sub <- array.sub[with(array.sub, order(gene, tissue, time)), ]
  
  # Check my gene, tissue and times all line up
  vars <- c("gene", "tissue", "time")
  print(sapply(vars, function(var) all(kallisto.sub[[var]] == array.sub[[var]])))
  
  # Now just cbind because the gene tissue and times line up
  array.signal <- array.sub[[array.var]]
  merge.long <- cbind(kallisto.sub, array.signal)
  return(merge.long)
}

GetWeightsFromVar <- function(dat, fit.noise){
  # from array signal, estimate the noise in the estimate from fit.noise
  # We can't have negative values, so take the smallest non-negative number
  # in the signal as a proxy
  
  # BEGIN: Set variance as weights
  M.var <- predict(fit.noise, dat$array.signal)
  # adjust M.var so all values less than 0 take on smallest non-negative number
  M.var.min <- min(M.var[which(M.var > 0)])
  # if M.var.min is infinity, probably RNA-Seq has no expression, default weights to 1
  if (is.infinite(M.var.min)){
    gene <- dat$gene[1]
    warning(paste("Gene:", gene, "...adjusting infinites to 1 for variance."))
    M.var.min <- 1
  }
  M.var[which(M.var < 0)] <- M.var.min
  weights <- 1 / M.var
  # END: Set variance as weights
  return(weights)
}

FitSaturationCurve <- function(dat, fit.noise, array.wide){
  if (sum(dat$tpm) == 0){
    # all zeros, then return linear slope with slope = infinity
    return(data.frame(model = "lm", 
                      a = NA, 
                      b = NA,
                      k = NA,
                      int = 1,
                      slope = Inf))
  }
  gene <- dat$gene[1]
  M.full <- array.wide[gene, ]
  # BEGIN: Constants
  # x.factor: make fit not so close to saturation points
  # x.factor = 1: best fit
  # x.factor > 1: force saturation points farhter from data points
  x.factor <- 1.2
  # bg.factor: make fit not so close to background
  bg.factor <- 0.8
  
  # lower and upper bounds: saturation model
  bmin <- 0
  kmin <- 0
  amin <- max(M.full) * x.factor
  bmax <- min(M.full) * bg.factor
  amax <- Inf
  kmax <- Inf
  # init vals: saturation model
  b0 <- 2^4
  k0 <- 2^6  # large to begin in linear regime
  a0 <- amin + 10  # expect it close to amin
  # lower and upper bounds: 
  slopemin <- 0
  intmin <- 0
  slopemax <- Inf
  intmax <- min(M.full) * bg.factor
  
  # Linear fit constraints
  slope0 <- 0.07  # may need adjusting
  int0 <- 10
  # END: Constants
  
  weights <- GetWeightsFromVar(dat, fit.noise)
  
  fit <- tryCatch({
    fit.saturation <- nls(data = dat,
                          formula = array.signal ~ b + (a * tpm) / (k + tpm),
                          algorithm = "port",
                          start=list(a=a0, 
                                     b=b0,
                                     k=k0),
                          lower=list(a=amin,
                                     b=bmin,
                                     k=kmin),
                          upper=list(a=amax,
                                     b=bmax,
                                     k=kmax),
                          weights=weights)
    
    #     fit.lm <- FitLmConstraint2(dat, weights,
    #                                int0, slope0,
    #                                intmin, slopemin,
    #                                intmax, slopemax)
    fit.lm <- c(int = NA, slope = NA)
    return(data.frame(model = "saturation", 
                      a = coefficients(fit.saturation)[["a"]], 
                      b = coefficients(fit.saturation)[["b"]],
                      k = coefficients(fit.saturation)[["k"]],
                      int = NA,
                      slope = NA))
  }, error = function(e){
    #     fit.saturation <- nls(data = dat,
    #                           formula = array.signal ~ b + (a * tpm) / (k + tpm),
    #                           algorithm = "port",
    #                           start=list(a=a0, 
    #                                      b=b0,
    #                                      k=k0),
    #                           lower=list(a=amin,
    #                                      b=bmin,
    #                                      k=kmin),
    #                           upper=list(a=amax,
    #                                      b=bmax,
    #                                      k=kmax),
    #                           weights=weights)
    fit.saturation <- c(a = NA, b = NA, k = NA)
    fit.lm <- FitLmConstraint2(dat, weights,
                               int0, slope0,
                               intmin, slopemin,
                               intmax, slopemax)
    return(data.frame(model = fit.lm[["method"]], 
                      a = NA, 
                      b = NA,
                      k = NA,
                      int = fit.lm[["int"]],
                      slope = fit.lm[["slope"]]))
  })
}

# Load matrices -----------------------------------------------------------

kallisto.path <- "data/kallisto/abundance.genecounts.matrix.txt"
array.path <- "data/array_exprs_colnames_fixed.best.probe.selected.txt"
kallisto.long <- LoadKallistoGene(kallisto.path)
array.long <- LoadArray(array.path, get.norm = TRUE, form = "long")
array.wide <- LoadArray(array.path, get.norm = TRUE, form = "wide")

# jgene <- "Elovl3"
# PlotGeneTpm(kallisto.long, jgene, log2.pseudocount=0.001, scale = 1)
# PlotGeneAcrossTissues(subset(dat.long, gene == jgene))


# # Density plot my TPM -----------------------------------------------------
# 
# scale.factor <- 100  # make kallisto range between 0 to 15 log2 scale when I do +1 pseudo count
# kallisto.long$tpm <- kallisto.long$tpm * scale.factor
# kallisto.sub <- subset(kallisto.long, tpm > 0)
# plot(density(log2(kallisto.sub$tpm + 1)))
# Get mean and variance for each gene across tissues ----------------------


# Get a model for my noise ------------------------------------------------

mean.var.gene.tiss <- array.long %>%
  group_by(gene, tissue) %>%
  summarise(Mean = mean(signal.norm),
            Var = var(signal.norm)) %>%
  .[order(.$Mean, decreasing = FALSE), ]

library(plyr)
n.per.bin <- 150
mean.var.gene.tiss$BinOrder <- factor(round_any(seq(1:nrow(mean.var.gene.tiss)), n.per.bin))
detach("package:plyr", unload=TRUE)
library(dplyr)

mean.var.meds <- mean.var.gene.tiss %>%
  group_by(BinOrder) %>%
  summarise(MedianExprs = median(Mean),
            MedianVar = median(Var))

# head(data.frame(mean.var.meds[order(mean.var.meds$MedianExprs), ]), n = 500)

fit.noise <- loess(formula = MedianVar ~ MedianExprs, 
                   data = mean.var.meds, 
                   control=loess.control(surface="direct"))

# Fit Loess Model
x <- seq(min(array.long$signal.norm), max(array.long$signal.norm), 100)  # plot full range
y <- predict(fit.noise, x)

# # Plot diagnostic plots
# pdf("plots/adjust_array_to_kallisto_diagnostics/loess_fit_diagnostics.pdf")
# plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
#      xlab="bin.mean", ylab="bin.var",
#      xlim=c(0, 3000), ylim=c(0, 50000))
# points(mean.var.meds$MedianExprs, mean.var.meds$MedianVar)
# 
# plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
#      xlab="bin.mean", ylab="bin.var")
# points(mean.var.meds$MedianExprs, mean.var.meds$MedianVar)
# 
# plot(x, y, col='red', lwd='2', type='l', main=paste0('Loess Fit. Bin size=', n.per.bin),
#      xlab="bin.mean", ylab="bin.var", log="xy")
# points(mean.var.meds$MedianExprs, mean.var.meds$MedianVar)
# dev.off()


# Merge RNA-Seq and array as a long data frame, separate columns ----------

merge.long <- MergeRnaseqArray(kallisto.long, array.long, array.var = "signal.norm")


fits.test <- subset(merge.long, gene %in% c(GetGeneList())) %>%
  group_by(gene) %>%
  do(FitSaturationCurve(., fit.noise, array.wide))
# 
test <- subset(merge.long, gene == "BC053393")
test.fit <- FitSaturationCurve(test, fit.noise, array.wide)

# # head(merge.long)
# fits.all <- merge.long %>%
#   group_by(gene) %>%
#   do(FitSaturationCurve(., fit.noise, array.wide))
# save(fits.all, file = "Robjs/adjust_array_to_kallisto.fits.all.Robj")


# Do some diagnostic plots ------------------------------------------------
# 
testgenes <- GetGeneList()
g <- testgenes[1]






# Fit saturation curve ----------------------------------------------------


