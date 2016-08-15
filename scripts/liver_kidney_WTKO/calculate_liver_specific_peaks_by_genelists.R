# 2016-08=08
# how many liver-specific peaks there are in different sets? Redo with liver kidney WT KO 

rm(list=ls())

jmc.cores <- 10

# jtiss <- "Heart"
args <- commandArgs(trailingOnly=TRUE)
# jtiss <- args[1]
# jtiss <- "Liver"
jmod <- args[1]  # Liver_SV129
jtiss <- args[2]  # Liver
n.trials <- as.numeric(args[3])
# jmod <- "Liver_SV129,Liver_BmalKO"
# jmod <- "Liver_SV129"

random.is.flat <- TRUE

print(paste("Tissue", jtiss))
print(paste("Module", jmod))
# n.trials <- 1000
if (is.na(n.trials)) stop("Trials not numeric")

start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(parallel)
library(hash)

# Functions ---------------------------------------------------------------

GetTissueSpecStats <- function(S.long, jgenes, jcutoff, jcutoff.low, distfilt, jlab, jtiss="Liver"){
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
    filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff.low)
  
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


distfilt <- 40000
jcutoff <- 3
jcutoff.low <- 0

# Load --------------------------------------------------------------------

load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
S.long <- subset(S.long, tissue %in% c("Liver", "Kidney"))
S.long$tissue <- factor(as.character(S.long$tissue), levels = c("Kidney", "Liver"))
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == "g=1001")

# Get liv genes -----------------------------------------------------------

jgenes.all <- as.character(fits.long.filt$gene)
jgenes <- as.character(subset(fits.long.filt, model == jmod)$gene)
# jgenes.mus <- as.character(subset(fits.best, model == "Mus")$gene)
jgenes.flat <- as.character(subset(fits.long.filt, model == "")$gene)
jgenes.flat.filt <- subset(fits.bytiss, gene %in% jgenes.flat) %>%
  group_by(gene) %>%
  filter(abs(diff(range(int))) <= 2)
jgenes.flat.filt <- unique(as.character(jgenes.flat.filt$gene))

# clock.controlled.genes <- as.character(subset(fits.all.long.wtkohog, model %in% c("WT,Liver", "WT;KO;Liver", "WT;Liver", "Liver"))$gene)
# jgenes.filt <- intersect(jgenes, clock.controlled.genes)

# DEBUG
# jgenes <- jgenes[1:20]
# jgenes.flat <- jgenes.flat[1:20]
# jgenes.flat.filt <- jgenes.flat.filt[1:20]

gene.lists <- list("jgenes"=jgenes, "jgenes.flat"=jgenes.flat, "jgenes.flat.filt"=jgenes.flat.filt)
if (jmc.cores == 1){
  df.out.lst <- lapply(gene.lists, function(gene.list){
    df.out <- GetTissueSpecStats(S.long, gene.list, jcutoff, jcutoff.low, distfilt, jtiss, jtiss = jtiss)
  })
} else {
  df.out.lst <- mclapply(gene.lists, function(gene.list){
    df.out <- GetTissueSpecStats(S.long, gene.list, jcutoff, jcutoff.low, distfilt, jtiss, jtiss = jtiss)
  }, mc.cores = jmc.cores)
}
df.out.lst <- do.call(rbind, df.out.lst)
df.out.lst$gene.type <- names(gene.lists)

print(df.out.lst)

# Add random genes --------------------------------------------------------

# set.seed(0)

gene.lists.rand <- list()
for (i in seq(n.trials)){
  if (random.is.flat){
    jgenes.rand <- sample(jgenes.flat.filt, length(jgenes))
  } else {
    jgenes.rand <- sample(jgenes.all, length(jgenes))
  }
  # jgenes.rand <- sample(jgenes.flat.filt, length(jgenes))
  gene.lists.rand[[i]] <- jgenes.rand
}
if (jmc.cores == 1){
  df.out.lst.rand <- lapply(gene.lists.rand, function(gene.list){
    df.out <- GetTissueSpecStats(S.long, gene.list, jcutoff, jcutoff.low, distfilt, "Random", jtiss = jtiss)
  })
} else {
  df.out.lst.rand <- mclapply(gene.lists.rand, function(gene.list){
    df.out <- GetTissueSpecStats(S.long, gene.list, jcutoff, jcutoff.low, distfilt, "Random", jtiss = jtiss)
  }, mc.cores = jmc.cores)
}
df.out.lst.rand <- do.call(rbind, df.out.lst.rand)
print(head(df.out.lst.rand))


# Calculate z-score -------------------------------------------------------

df.out.lst.merged <- rbind(as.data.frame(df.out.lst), as.data.frame(df.out.lst.rand))
robjdir <- "/home/yeung/projects/tissue-specificity/Robjs/liver_kidney_atger_nestle/tissue_specific_peaks"
outf <- file.path(robjdir, paste0("n.tiss.spec.df.out.lst.rand.", n.trials, ".tissue.", jtiss, ".module.", jmod, ".random.flat.", random.is.flat, ".Robj"))
# save(df.out.lst.merged, file = "Robjs/n.liv.spec.df.out.lst.rand.Robj")
save(df.out.lst.merged, file = outf)
print(Sys.time() - start)

