# plot_exprs_from_genelist.R
# March 3 2015
library(ggplot2)
library(plyr)
library(foreach)
library(doParallel)
setwd("/home/yeung/projects/tissue-specificity")

# Functions ---------------------------------------------------------------

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

# Set cluster -------------------------------------------------------------

ncores <- detectCores()
registerDoParallel(cores = ncores)


# Load data ---------------------------------------------------------------


scripts.dir <- "scripts"
funcs.dir <- "functions"
source(file.path(scripts.dir, funcs.dir, "LoadArrayRnaSeq.R"))
source(file.path(scripts.dir, funcs.dir, "FitRhythmic.R"))
source(file.path(scripts.dir, funcs.dir, "PlotGeneAcrossTissues.R"))
# source(file.path(scripts.dir, funcs.dir, "RemoveExtension.R"))
# source(file.path(scripts.dir, funcs.dir, "ReadListToVector.R"))

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
tissue_specific_genes <- subset(tissue_specific_genes, IsTissueSpecific == TRUE)


# Doing SVD on this matrix. What happens? ---------------------------------


