# 2016-02-19
# do penalized LDA using DHS on gene body
# redo using a new set of sitecounts matrix derived from get_tissue_spec_peaks.R

rm(list=ls())

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/dat.fit.Robj", v=T)
source("scripts/functions/PlotGeneAcrossTissues.R")

load("Robjs/N.long.sum.bytiss.all_genes.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)



library(reshape2)
library(penalizedLDA)



# Functions ---------------------------------------------------------------

RunPenalizedLDA <- function(N.long.sum.bytiss, fits.best, jmodels, fg.tiss, bg.tiss, bg.genes.type, sitecount.name, minamp, n.bg.genes = "all"){
  
  if (jmodels != ""){
    gene.list <- as.character(subset(fits.best, model %in% jmodels & amp.avg > minamp)$gene)
  } else {
    # ignore amplitude requirement if model is flat
    gene.list <- as.character(subset(fits.best, model %in% jmodels)$gene)
  }
  
  # define flat genes
  if (bg.genes.type == "Flat"){
    gene.list.bg <- as.character(subset(fits.best, model == "")$gene)
    #  gene.list.bg <- unique(as.character(subset(dat.fit, tissue == "Liver" & gene %in% gene.list.bg & int.rnaseq > 9)$gene))
  } else if (bg.genes.type == "Rhythmic"){
    gene.list.bg <- as.character(subset(fits.best, n.rhyth >= 8 & amp.avg > minamp)$gene)
  } else if (bg.genes.type == "same"){
    gene.list.bg <- gene.list  # use same genes
  } else{
    warning("Unexpected bg.genes.type")
  }
  
  set.seed(0)
  if (n.bg.genes == "all"){
    # take all bg genes
    print("No subsampling for bg genes")
  } else if (n.bg.genes == "same"){
    n.genes <- length(gene.list)
    gene.list.bg <- sample(gene.list.bg, size = n.genes)
  } else {
    # assume a numeric
    gene.list.bg <- sample(gene.list.bg, size = n.bg.genes)
  }
  
  
  N.sub <- subset(N.long.sum.bytiss, gene %in% gene.list & tissue %in% fg.tiss)
  N.bkgrd.all <- subset(N.long.sum.bytiss, gene %in% gene.list.bg & tissue %in% bg.tiss)
  # N.bkgrd.all <- subset(N.long.sum.bytiss, gene %in% gene.list.bg & tissue %in% fg.tiss)
  
  # sample to have same number as N.sub
  (n.frgrd <- nrow(N.sub))
  gene.list.N <- unique(N.sub$gene)
  
  # Run LDA -----------------------------------------------------------------
  
  # make a new label gene;tissue which will be my rowname
  N.sub$genetiss <- paste(N.sub$gene, N.sub$tissue, sep = ";")
  N.bkgrd.all$genetiss <- paste(N.bkgrd.all$gene, N.bkgrd.all$tissue, sep = ";")
  
  N.sub.mat <- dcast(data = N.sub, formula = genetiss ~ motif, fill = 0, value.var = sitecount.name)
  N.bkgrd.all.mat <- dcast(N.bkgrd.all, genetiss ~ motif, fill = 0, value.var = sitecount.name)
  
  common.cnames <- intersect(colnames(N.sub.mat), colnames(N.bkgrd.all.mat))
  # ignore non-common cnames
  N.sub.mat <- N.sub.mat[, common.cnames]
  N.bkgrd.all.mat <- N.bkgrd.all.mat[, common.cnames]
  
  labs <- c(rep(1, nrow(N.sub.mat)), rep(2, nrow(N.bkgrd.all.mat)))
  
  N.mat.merged <- rbind(N.sub.mat, N.bkgrd.all.mat)
  rownames(N.mat.merged) <- N.mat.merged$genetiss; N.mat.merged$genetiss <- NULL
  
  out <- PenalizedLDA(as.matrix(N.mat.merged), labs, lambda = 0.1, K = 1, standardized = FALSE)
  
  
  out.df <- data.frame(xproj = out$xproj, lab = out$y)
  boxplot(xproj ~ lab, data = out.df)
  
  discrim.filt <- sort(out$discrim[which(out$discrim != 0)], decreasing = FALSE)
  cnames <- colnames(N.mat.merged)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
  plot(seq(length(discrim.filt)), discrim.filt, 
       main=paste(paste(fg.tiss, collapse = ","), "vs", paste(bg.tiss, collapse=","), sitecount.name, "\n", 
                  "ngenes fg:", length(gene.list), "ngenes bg:", length(gene.list.bg), "\n",
                  "models", paste0(jmodels, collapse=","), "\n",
                  "bg genes type:", bg.genes.type))
  text(seq(length(discrim.filt)), discrim.filt, labels = cnames)
  
  print(cnames)
  return(out)
}


# Load --------------------------------------------------------------------

minamp <- 0.25

# define parameters
# use.cross <- FALSE
sitecount.name <- "sitecount.max"
# bg.genes.type <- "Rhythmic"  # "Flat", "Rhythmic", "same" 
bg.genes.type <- "Flat"  # "Flat", "Rhythmic", "same" 

# define rhythmic model
# jmodels <- c("Kidney")

# Kidney Liver models
jgrep <- "^Kidney;Liver$|^Kidney,Liver$"
jmodels <- unique(fits.best[grepl(jgrep, fits.best$model), ]$model)

# define foreground and background tissues
all.tiss <- c('Cere', 'Heart', 'Kidney', 'Liver', 'Lung', 'Mus')
fg.tiss <- c("Liver", "Kidney")
# fg.tiss <- c("Liver")
# fg.tiss <- c("Kidney")
# fg.tiss <- c("Heart")
# bg.tiss <- all.tiss[which(! all.tiss %in% fg.tiss)]
bg.tiss <- all.tiss[which(! all.tiss %in% fg.tiss & all.tiss != "Kidney")]
# bg.tiss <- c("Kidney")
# bg.tiss <- c("Liver")
print(fg.tiss)
print(bg.tiss)

RunPenalizedLDA(N.long.sum.bytiss, fits.best, jmodels, fg.tiss, fg.tiss, bg.genes.type, sitecount.name, minamp = 0.25, n.bg = 160)


# 
# if (jmodels != ""){
#   gene.list <- as.character(subset(fits.best, model %in% jmodels & amp.avg > minamp)$gene)
# } else {
#   # ignore amplitude requirement if model is flat
#   gene.list <- as.character(subset(fits.best, model %in% jmodels)$gene)
# }
# 
# # define flat genes
# if (bg.genes.type == "Flat"){
#   gene.list.bg <- as.character(subset(fits.best, model == "")$gene); gene.list.bg <- unique(as.character(subset(dat.fit, tissue == "Liver" & gene %in% gene.list.bg & int.rnaseq > 9)$gene))
# } else if (bg.genes.type == "Rhythmic"){
#   gene.list.bg <- as.character(subset(fits.best, n.rhyth >= 8 & amp.avg > minamp)$gene)
# } else if (bg.genes.type == "same"){
#   gene.list.bg <- gene.list  # use same genes
# } else{
#   warning("Unexpected bg.genes.type")
# }
# 
# N.sub <- subset(N.long.sum.bytiss, gene %in% gene.list & tissue %in% fg.tiss)
# N.bkgrd.all <- subset(N.long.sum.bytiss, gene %in% gene.list.bg & tissue %in% bg.tiss)
# # N.bkgrd.all <- subset(N.long.sum.bytiss, gene %in% gene.list.bg & tissue %in% fg.tiss)
# 
# # sample to have same number as N.sub
# (n.frgrd <- nrow(N.sub))
# gene.list.N <- unique(N.sub$gene)
# 
# # Run LDA -----------------------------------------------------------------
# 
# # make a new label gene;tissue which will be my rowname
# N.sub$genetiss <- paste(N.sub$gene, N.sub$tissue, sep = ";")
# N.bkgrd.all$genetiss <- paste(N.bkgrd.all$gene, N.bkgrd.all$tissue, sep = ";")
# 
# 
# N.sub.mat <- dcast(data = N.sub, formula = genetiss ~ motif, fill = 0, value.var = sitecount.name)
# N.bkgrd.all.mat <- dcast(N.bkgrd.all, genetiss ~ motif, fill = 0, value.var = sitecount.name)
# 
# common.cnames <- intersect(colnames(N.sub.mat), colnames(N.bkgrd.all.mat))
# # ignore non-common cnames
# N.sub.mat <- N.sub.mat[, common.cnames]
# N.bkgrd.all.mat <- N.bkgrd.all.mat[, common.cnames]
# 
# labs <- c(rep(1, nrow(N.sub.mat)), rep(2, nrow(N.bkgrd.all.mat)))
# 
# N.mat.merged <- rbind(N.sub.mat, N.bkgrd.all.mat)
# rownames(N.mat.merged) <- N.mat.merged$genetiss; N.mat.merged$genetiss <- NULL
# 
# out <- PenalizedLDA(as.matrix(N.mat.merged), labs, lambda = 0.1, K = 1, standardized = FALSE)
# 
# print(out)
# 
# out.df <- data.frame(xproj = out$xproj, lab = out$y)
# boxplot(xproj ~ lab, data = out.df)
# 
# discrim.filt <- sort(out$discrim[which(out$discrim != 0)], decreasing = FALSE)
# cnames <- colnames(N.mat.merged)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
# plot(seq(length(discrim.filt)), discrim.filt, 
#      main=paste(paste(fg.tiss, collapse = ","), "vs", paste(bg.tiss, collapse=","), sitecount.name, "\n", 
#                 "ngenes fg:", length(gene.list), "ngenes bg:", length(gene.list.bg)))
# text(seq(length(discrim.filt)), discrim.filt, labels = cnames)
# 
# print(cnames[1:20])
# 
# # 
# # plot(out$xproj, out$y)
# # 
# # plot separation using a single dimension
# # jmotif <- "RXRG_dimer.p3"
# # jmotif <- "ONECUT1,2.p2"
# # jmotif <- "HNF1A.p2"
# # jmotif <- "FOXA2.p3"
# # jmotif <- "MAFB.p2"
# # 
# # jmotif <- "bHLH_family.p2"
# # out.df.motif <- data.frame(sitecount = N.mat.merged[, jmotif], lab = labs)
# # boxplot(sitecount ~ lab, data = out.df.motif, main = jmotif)
# # plot(labs, N.mat.merged[, jmotif])
# 
# 
# # Where are the RXRG_dimer motifs located? --------------------------------
# 
# 
