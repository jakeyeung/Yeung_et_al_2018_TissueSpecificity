# 2016-04-22
# Subset plda, also ask how many liver-specific peaks there are in different sets?

rm(list=ls())

start <- Sys.time()
setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)
library(parallel)
library(hash)

# Functions ---------------------------------------------------------------

GetLiverSpecStats <- function(S.long, jgenes, jcutoff, distfilt, jlab){
  # Mean number of liver-specific DHS per gene
  # Total genes with liver-specific DHS
  # Total number of liver-specific DHS
  # Total genes
  # Fractio nof genes with liver-specific DHS
  
  # get peaks near gene
  S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
  jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
  
  # Identify tissue-specific peaks ------------------------------------------
  
  rhyth.tiss <- c("Liver")
  
  # take peaks with Liver signal greater than cutoff
  jtiss <- levels(S.sub$tissue)
  tiss.i <- which(jtiss %in% rhyth.tiss)
  others.i <- which(!jtiss %in% rhyth.tiss)
  
  S.sub.liverpeaks <- S.sub %>%
    group_by(peak) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
    # filter(mean(zscore[tiss.i]) > jcutoff & mean(zscore[others.i]) < jcutoff)
    filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff)
  
  # Find how many peaks -----------------------------------------------------
  
  # label peaks as liver-specific or not
  
  liv.spec.peaks <- hash(as.character(S.sub.liverpeaks$peak), TRUE)
  
  S.sub$is.liv.spec <- sapply(as.character(S.sub$peak), function(p){
    is.liv.spec <- liv.spec.peaks[[p]]
    if (is.null(is.liv.spec)){
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  S.sub.n.livspec <- subset(S.sub, tissue == "Liver") %>%
    group_by(tissue, gene) %>%
    summarise(n.liv.spec = length(which(is.liv.spec == TRUE))) %>%
    mutate(gene.type = jlab)
  
  
  # How many genes have liver-specific genes? -------------------------------
  
  n.liv.spec.sum <- S.sub.n.livspec %>%
    group_by(gene.type) %>%
    summarise(mean.liv.spec.per.gene = sum(n.liv.spec) / length(n.liv.spec), 
              total.liv.spec = sum(n.liv.spec), 
              total.genes = length(n.liv.spec), 
              genes.with.liv.spec = length(which(n.liv.spec > 0)), 
              frac.n.spec.by.gene = genes.with.liv.spec / total.genes)
  return(n.liv.spec.sum)
}


# Constnts ----------------------------------------------------------------


distfilt <- 5000
jcutoff <- 1.5
# Load --------------------------------------------------------------------

load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)
load("Robjs/N.long.livertwflat.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
load("Robjs/wtko.dat.wtko.hog.commongenes.Robj", v=T)
load("Robjs/wtko.fits.all.long.commongenes.Robj", v=T)
load("Robjs/dat.fit.Robj", v=T)

# Get liv genes -----------------------------------------------------------

jgenes.all <- as.character(fits.best$gene)
jgenes <- as.character(subset(fits.best, model == "Liver")$gene)
jgenes.flat <- as.character(subset(fits.best, model == "")$gene)
jgenes.flat.filt <- subset(dat.fit, gene %in% jgenes.flat) %>%
  group_by(gene) %>%
  filter(abs(diff(range(int.rnaseq))) <= 2)
jgenes.flat.filt <- unique(as.character(jgenes.flat.filt$gene))

clock.controlled.genes <- as.character(subset(fits.all.long.wtkohog, model %in% c("WT,Liver", "WT;KO;Liver", "WT;Liver", "Liver"))$gene)
jgenes.filt <- intersect(jgenes, clock.controlled.genes)

gene.lists <- list("jgenes"=jgenes, "jgenes.filt"=jgenes.filt, "jgenes.flat"=jgenes.flat, "jgenes.flat.filt"=jgenes.flat.filt)
df.out.lst <- mclapply(gene.lists, function(gene.list){
  df.out <- GetLiverSpecStats(S.long, gene.list, jcutoff, distfilt, "TODO")
}, mc.cores = 4)
df.out.lst <- do.call(rbind, df.out.lst)
df.out.lst$gene.type <- names(gene.lists)


# Add random genes --------------------------------------------------------

set.seed(0)

gene.lists.rand <- list()
n.trials <- 10000
for (i in seq(n.trials)){
  # jgenes.rand <- sample(jgenes.all, length(jgenes))
  jgenes.rand <- sample(jgenes.flat.filt, length(jgenes))
  gene.lists.rand[[i]] <- jgenes.rand
}
df.out.lst.rand <- mclapply(gene.lists.rand, function(gene.list){
  df.out <- GetLiverSpecStats(S.long, gene.list, jcutoff, distfilt, "Random")
}, mc.cores = 10)
print(head(df.out.lst.rand))
df.out.lst.rand <- do.call(rbind, df.out.lst.rand)


# Calculate z-score -------------------------------------------------------

df.out.lst.merged <- rbind(df.out.lst, df.out.lst.rand)
outf <- paste0("Robjs/n.liv.spec.df.out.lst.rand.similarlyexpressed.", n.trials, ".Robj")
# save(df.out.lst.merged, file = "Robjs/n.liv.spec.df.out.lst.rand.Robj")
save(df.out.lst.merged, file = outf)

print(Sys.time() - start)


# Downstream --------------------------------------------------------------

do.downstream <- FALSE
if (do.downstream){
  load("Robjs/n.liv.spec.df.out.lst.rand.Robj", v=T)
  ggplot(subset(df.out.lst.merged, gene.type == "Random"), aes(x = frac.n.spec.by.gene)) + geom_histogram() + geom_vline(aes(xintercept = 0.43)) + ggtitle("Fraction of genes with at least one Liver-specific DHS (random)")
}

# ggplot(df.out.lst.rand, aes(x = frac.n.spec.by.gene)) + geom_histogram(bins = 50)
# ggplot(df.out.lst.rand, aes(x = frac.n.spec.by.gene)) + geom_density()


# 
# df.out <- GetLiverSpecStats(S.long, gene.list, jcutoff, distfilt, "RhythLiver")
# 
# # S.sub.liv <- subset(S.long, gene %in% liv.genes)
# # S.sub.flat <- subset(S.long, gene %in% flat.genes)
# 
# # get peaks near gene
# S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
# S.sub.flat <- subset(S.long, gene %in% jgenes.flat & dist < distfilt)
# jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
# jpeaks.flat <- as.character(unique(S.sub.flat$peak))
# 
# print(paste("number of peaks in liver-specific rhythmic genes", length(jpeaks)))
# print(paste("number of peaks in flat genes", length(jpeaks.flat)))
# 
# N.sub <- subset(N.long.filt, gene %in% jgenes & peak %in% jpeaks)  # should not have any missing peaks
# N.sub.flat <- subset(N.long.filt, gene %in% jgenes.flat & peak %in% jpeaks.flat)
# 
# 
# # Identify tissue-specific peaks ------------------------------------------
# 
# 
# rhyth.tiss <- c("Liver")
# 
# # take peaks with Liver signal greater than cutoff
# jtiss <- levels(S.sub$tissue)
# tiss.i <- which(jtiss %in% rhyth.tiss)
# others.i <- which(!jtiss %in% rhyth.tiss)
# 
# S.sub.liverpeaks <- S.sub %>%
#   group_by(peak) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
#   # filter(mean(zscore[tiss.i]) > jcutoff & mean(zscore[others.i]) < jcutoff)
#   filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff)
# 
# jtiss.flat <- levels(S.sub.flat$tissue)
# if (identical(jtiss, jtiss.flat) == FALSE){
#   print("This shouldnt be necessary")
#   tiss.i <- which(jtiss %in% rhyth.tiss)
#   others.i <- which(jtiss %in% rhyth.tiss)
# }
# 
# S.sub.flat.liverpeaks <- S.sub.flat %>%
#   group_by(peak) %>%
#   # filter(mean(zscore[tiss.i]) > jcutoff & mean(zscore[others.i]) < jcutoff)
#   filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff)
# 
# 
# # Find how many peaks -----------------------------------------------------
# 
# # label peaks as liver-specific or not
# 
# liv.spec.peaks <- hash(as.character(S.sub.liverpeaks$peak), TRUE)
# liv.spec.peaks.flat <- hash(as.character(S.sub.flat.liverpeaks$peak), TRUE)
# 
# S.sub$is.liv.spec <- sapply(as.character(S.sub$peak), function(p){
#   is.liv.spec <- liv.spec.peaks[[p]]
#   if (is.null(is.liv.spec)){
#     return(FALSE)
#   } else {
#     return(TRUE)
#   }
# })
# S.sub.n.livspec <- subset(S.sub, tissue == "Liver") %>%
#   group_by(tissue, gene) %>%
#   summarise(n.liv.spec = length(which(is.liv.spec == TRUE))) %>%
#   mutate(gene.type = "RhythLiver")
# 
# S.sub.flat$is.liv.spec <- sapply(as.character(S.sub.flat$peak), function(p){
#   is.liv.spec <- liv.spec.peaks.flat[[p]]
#   if (is.null(is.liv.spec)){
#     return(FALSE)
#   } else {
#     return(TRUE)
#   }
# })
# 
# S.sub.flat.n.livspec <- subset(S.sub.flat, tissue == "Liver") %>%
#   group_by(tissue, gene) %>%
#   summarise(n.liv.spec = length(which(is.liv.spec == TRUE))) %>%
#   mutate(gene.type = "Flat")
# 
# 
# # How many genes have liver-specific genes? -------------------------------
# 
# S.sub.n.livspec.merged <- rbind(S.sub.n.livspec, S.sub.flat.n.livspec)
# 
# n.liv.spec.sum <- S.sub.n.livspec.merged %>%
#   group_by(gene.type) %>%
#   summarise(mean.liv.spec.per.gene = sum(n.liv.spec) / length(n.liv.spec), 
#             total.liv.spec = sum(n.liv.spec), 
#             total.genes = length(n.liv.spec), 
#             genes.with.liv.spec = length(which(n.liv.spec > 0)), 
#             frac.n.spec.by.gene = genes.with.liv.spec / total.genes)
# 
# # visualize by distribution
# ggplot(S.sub.n.livspec.merged, aes(x = n.liv.spec)) + geom_histogram() + facet_wrap(~gene.type)

