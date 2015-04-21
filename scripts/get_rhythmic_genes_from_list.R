# plot_exprs_from_genelist.R
# March 3 2015
library(ggplot2)
library(plyr)
library(foreach)
library(doParallel)
setwd("/home/yeung/projects/tissue-specificity")


# Constants ---------------------------------------------------------------

T <- 24  # hours
omega <- 2 * pi / T

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/ConvertLongToWide.R")
source("scripts/functions/PlotActivitiesFunctions.R")

Is.TissueSpecific <- function(pval, amp, min.pval = 1e-5, max.pval = 0.05, max.amp = 0.5, min.amp = 0.1){
  # Check that min pval and max amplitude passes cutoff
  # And max pval and min amplitude passes cutoff
  
  # first two checks: if there exists a rhythmic gene
  # last two checks: if there exists a non-rhythmic gene.
  if (min(pval) < min.pval & max(amp) > max.amp & 
        max(pval) > max.pval & min(amp) < min.amp){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

Is.RhythmicAcrossTissues <- function(pval, amp, cutoff.pval = 1e-5, cutoff.amp = 0.5){
  if (max(pval) < cutoff.pval & min(amp) > cutoff.amp){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Set cluster -------------------------------------------------------------

ncores <- detectCores()
registerDoParallel(cores = ncores)


# Load data ---------------------------------------------------------------


scripts.dir <- "scripts"
funcs.dir <- "functions"

genes <- ReadListToVector("data/gene_lists/genes_associated_with_multiple_promoters.txt")
genes <- ReadListToVector("data/gene_lists/genes_with_multiple_promoters_aftermm10liftOver.txt")
dat <- LoadArrayRnaSeq()
dat.sub <- subset(dat, gene %in% genes)

# Fit Rhythmic ------------------------------------------------------------


# Fit rhythmic
sprintf("Running with %s cores.", ncores)
t.start <- Sys.time()
dat.sub.split <- split(dat.sub, dat.sub$tissue)
# remove WFAT
dat.sub.split$WFAT <- NULL
dat.sub.fit <- lapply(dat.sub.split, function(df){
  ddply(df, .(gene), FitRhythmic, .parallel=TRUE)
})
print(t.start - Sys.time())

# Create long dataframe
dat.sub.fit.long <- do.call(rbind, dat.sub.fit)

# Summarize each gene by average and SD of pval
cutoff.pval <- 1e-5
cutoff.amp <- 0.5

# Filter out neural tissues
# Filter out WFAT because it's a weird one.
dat.sub.fit.long.filt <- subset(dat.sub.fit.long, ! tissue %in% c("BS", "Cere", "Hypo", "WFAT"))


# Filter for genes containing significant pvals AND non significant across tissues
tissue_specific_genes <- ddply(dat.sub.fit.long.filt, .(gene), summarise,
                               IsTissueSpecific = Is.TissueSpecific(pval, amp))

tissue_specific_genes <- subset(tissue_specific_genes, IsTissueSpecific == TRUE)$gene

rhythmic_genes <- ddply(dat.sub.fit.long.filt, .(gene), summarise,
                        IsRhythmic = Is.RhythmicAcrossTissues(pval, amp))


# Doing SVD on this matrix. What happens? ---------------------------------

dat.sub.filt <- subset(dat.sub, gene %in% tissue_specific_genes)
dat.sub.split <- split(dat.sub.filt, dat.sub.filt$tissue)

dat.split.proj <- lapply(dat.sub.split, function(x){ 
  ddply(x, .(gene), ProjectToFrequency2, omega = omega, add.tissue = TRUE, .parallel = TRUE)
})
dat.proj <- do.call(rbind, dat.split.proj)



# Make wide format --------------------------------------------------------

dat.proj.wide <- ConvertLongToWide(dat.proj)



# Perform SVD -------------------------------------------------------------

tissues <- colnames(dat.proj.wide)
genes <- rownames(dat.proj.wide)

s <- svd(dat.proj.wide)
rownames(s$u) <- genes
rownames(s$v) <- tissues

plot(s$d^2 / sum(s$d ^ 2), type = 'o')

for (comp in seq(11)){
  comp <- 4
  PlotEigengene(s, comp)
  PlotEigensamp(s, comp)
}

