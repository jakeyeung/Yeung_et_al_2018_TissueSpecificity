# 2016-04-22
# Subset plda, also ask how many liver-specific peaks there are in different sets?

rm(list=ls())

jmc.cores <- 1

# jtiss <- "Heart"
args <- commandArgs(trailingOnly=TRUE)
jtiss <- args[1]
print(paste("Tissue", jtiss))
n.trials <- as.numeric(args[2])
if (is.na(n.trials)) stop("Trials not numeric")

start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(parallel)
library(hash)

# Functions ---------------------------------------------------------------

GetTissueSpecStats <- function(S.long, jgenes, jcutoff, distfilt, jlab, jtiss="Liver"){
  # Mean number of liver-specific DHS per gene
  # Total genes with liver-specific DHS
  # Total number of liver-specific DHS
  # Total genes
  # Fractio nof genes with liver-specific DHS
  
  # get peaks near gene
  S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
  jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
  
  # Identify tissue-specific peaks ------------------------------------------
  
  rhyth.tiss <- c(jtiss)
  
  # take peaks with Liver signal greater than cutoff
  tiss.all <- levels(S.sub$tissue)
  tiss.i <- which(tiss.all %in% rhyth.tiss)
  others.i <- which(!tiss.all %in% rhyth.tiss)
  if (length(tiss.i) == 0) stop(paste("Unknown tissue"))
  
  S.sub.tissuepeaks <- S.sub %>%
    group_by(peak) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
    # filter(mean(zscore[tiss.i]) > jcutoff & mean(zscore[others.i]) < jcutoff)
    filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff)
  
  # Find how many peaks -----------------------------------------------------
  
  # label peaks as tissue-specific or not
  jpeaks <- as.character(S.sub.tissuepeaks$peak)
  if (length(jpeaks > 1)){
    tiss.spec.peaks <- hash(as.character(S.sub.tissuepeaks$peak), TRUE)
  } else {
    tiss.spec.peaks <- hash("nopeaks", TRUE)
  }
  
  
  S.sub$is.tiss.spec <- sapply(as.character(S.sub$peak), function(p){
    is.tiss.spec <- tiss.spec.peaks[[p]]
    if (is.null(is.tiss.spec)){
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  S.sub.n.livspec <- subset(S.sub, tissue == "Liver") %>%
    group_by(tissue, gene) %>%
    summarise(n.tiss.spec = length(which(is.tiss.spec == TRUE))) %>%
    mutate(gene.type = jlab)
  
  
  # How many genes have tissue-specific genes? -------------------------------
  
  n.tiss.spec.sum <- S.sub.n.livspec %>%
    group_by(gene.type) %>%
    summarise(mean.tiss.spec.per.gene = sum(n.tiss.spec) / length(n.tiss.spec), 
              total.tiss.spec = sum(n.tiss.spec), 
              total.genes = length(n.tiss.spec), 
              genes.with.tiss.spec = length(which(n.tiss.spec > 0)), 
              frac.n.spec.by.gene = genes.with.tiss.spec / total.genes)
  return(n.tiss.spec.sum)
}


# Constnts ----------------------------------------------------------------


distfilt <- 5000
jcutoff <- 1.5
# Load --------------------------------------------------------------------

load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
# load("Robjs/N.long.livertwflat.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
load("Robjs/wtko.dat.wtko.hog.commongenes.Robj", v=T)
load("Robjs/wtko.fits.all.long.commongenes.Robj", v=T)
load("Robjs/dat.fit.Robj", v=T)

# Get liv genes -----------------------------------------------------------

jgenes.all <- as.character(fits.best$gene)
jgenes <- as.character(subset(fits.best, model == jtiss)$gene)
# jgenes.mus <- as.character(subset(fits.best, model == "Mus")$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)
jgenes.flat.filt <- subset(dat.fit, gene %in% jgenes.flat) %>%
  group_by(gene) %>%
  filter(abs(diff(range(int.rnaseq))) <= 2)
jgenes.flat.filt <- unique(as.character(jgenes.flat.filt$gene))

# clock.controlled.genes <- as.character(subset(fits.all.long.wtkohog, model %in% c("WT,Liver", "WT;KO;Liver", "WT;Liver", "Liver"))$gene)
# jgenes.filt <- intersect(jgenes, clock.controlled.genes)

gene.lists <- list("jgenes"=jgenes, "jgenes.flat"=jgenes.flat, "jgenes.flat.filt"=jgenes.flat.filt)
if (jmc.cores == 1){
  df.out.lst <- lapply(gene.lists, function(gene.list){
    df.out <- GetTissueSpecStats(S.long, gene.list, jcutoff, distfilt, jtiss, jtiss = jtiss)
  })
} else {
  df.out.lst <- mclapply(gene.lists, function(gene.list){
    df.out <- GetTissueSpecStats(S.long, gene.list, jcutoff, distfilt, jtiss, jtiss = jtiss)
  }, mc.cores = jmc.cores)
}
df.out.lst <- do.call(rbind, df.out.lst)
df.out.lst$gene.type <- names(gene.lists)

print(df.out.lst)

# Add random genes --------------------------------------------------------

set.seed(0)

gene.lists.rand <- list()
for (i in seq(n.trials)){
  jgenes.rand <- sample(jgenes.all, length(jgenes))
  # jgenes.rand <- sample(jgenes.flat.filt, length(jgenes))
  gene.lists.rand[[i]] <- jgenes.rand
}
if (jmc.cores == 1){
  df.out.lst.rand <- lapply(gene.lists.rand, function(gene.list){
    df.out <- GetTissueSpecStats(S.long, gene.list, jcutoff, distfilt, "Random", jtiss = jtiss)
  })
} else {
  df.out.lst.rand <- mclapply(gene.lists.rand, function(gene.list){
    df.out <- GetTissueSpecStats(S.long, gene.list, jcutoff, distfilt, "Random", jtiss = jtiss)
  }, mc.cores = jmc.cores)
}
df.out.lst.rand <- do.call(rbind, df.out.lst.rand)
print(head(df.out.lst.rand))


# Calculate z-score -------------------------------------------------------

df.out.lst.merged <- rbind(df.out.lst, df.out.lst.rand)
outf <- paste0("Robjs/n.tiss.spec.df.out.lst.rand.", n.trials, ".tissue.", jtiss, ".Robj")
# save(df.out.lst.merged, file = "Robjs/n.liv.spec.df.out.lst.rand.Robj")
save(df.out.lst.merged, file = outf)
print(Sys.time() - start)

