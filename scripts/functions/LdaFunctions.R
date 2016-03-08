library(penalizedLDA)
library(wordcloud)

DoLdaPromoters <- function(fg.mat, bg.mat){
  # remove columns with 0 within-class SD
  fg.mat <- fg.mat[, apply(fg.mat, 2, sd) > 0]
  bg.mat <- bg.mat[, apply(bg.mat, 2, sd) > 0]
  #   flat.mat <- flat.mat[, apply(flat.mat, 2, sd) > 0]
  #   rhyth.mat <- rhyth.mat[, apply(flat.mat, 2, sd) > 0]
  
  fg.labs <- rep(1, nrow(fg.mat))
  bg.labs <- rep(2, nrow(bg.mat))
  
  common.cols <- intersect(colnames(fg.mat), colnames(bg.mat))
  
  fg.mat <- fg.mat[, common.cols]
  bg.mat <- bg.mat[, common.cols]
  
  merged.mat <- rbind(fg.mat, bg.mat)
  merged.labs <- c(fg.labs, bg.labs)
  
  out <- PenalizedLDA(x = merged.mat, y = merged.labs, lambda = 0.11, K = 1)
  return(out)
}

LongToMat <- function(N.long, jvar = "sitecount"){
  mat <- dcast(N.long, formula = gene.uniq ~ motif, value.var = jvar)
  rownames(mat) <- mat$gene.uniq; mat$gene.uniq <- NULL
  return(mat)
}

RunPenalizedLDA <- function(N.long.sum.bytiss, fits.best, jmodels, fg.tiss, bg.tiss, bg.genes.type, sitecount.name, minamp, n.bg.genes = "all", jlambda=0.1){
  # from sitecounts_analysis_dhs_on_gene_body.redo.R
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
    # print("No subsampling for bg genes")
    
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
  
  out <- PenalizedLDA(as.matrix(N.mat.merged), labs, lambda = jlambda, K = 1, standardized = FALSE)  
  
#   out.df <- data.frame(xproj = out$xproj, lab = out$y)
#   boxplot(xproj ~ lab, data = out.df)
  
#   discrim.filt <- sort(out$discrim[which(out$discrim != 0)], decreasing = FALSE)
#   cnames <- colnames(N.mat.merged)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
#   plot(seq(length(discrim.filt)), discrim.filt, 
#        main=paste(paste(fg.tiss, collapse = ","), "vs", paste(bg.tiss, collapse=","), sitecount.name, "\n", 
#                   "ngenes fg:", length(gene.list), "ngenes bg:", length(gene.list.bg), "\n",
#                   "models", paste0(jmodels, collapse=","), "\n",
#                   "bg genes type:", bg.genes.type))
#   text(seq(length(discrim.filt)), discrim.filt, labels = cnames)
  out$fg.genes <- gene.list
  return(out)
}

PenalizedLdaLong <- function(fits.best, N.long, jmodel, jmodel.bg, jlambda, K, jvalue.var = "sitecount"){
  # fits.best: from Nconds
  # N.long: sitecounts from promoters
  # jmodel: get genes belonging to jmodel (can be vector) acts as foreground
  # jmodel.bg: background model for comparing against jmodel
  
  genes.fg <- subset(fits.best, model %in% jmodel)$gene
  genes.bg <- subset(fits.best, model %in% jmodel.bg)$gene
  
  N.fg <- subset(N.long, gene %in% genes.fg)
  N.bg <- subset(N.long, gene %in% genes.bg)
  
  N.mat.fg <- dcast(data = subset(N.long, gene %in% genes.fg), 
                    formula = gene.uniq ~ motif,
                    value.var = jvalue.var,
                    fill = 0)
  N.mat.bg <- dcast(data = subset(N.long, gene %in% genes.bg), 
                    formula = gene.uniq ~ motif,
                    value.var = jvalue.var,
                    fill = 0)
  
  N.mat.fgbg <- rbind(N.mat.fg, N.mat.bg)
  # add labels starting from 1...  to correspond to fg or bg
  labs <- c(rep(1, nrow(N.mat.fg)), rep(2, nrow(N.mat.bg)))
  
  # remove first column name move it to rownames
  rownames(N.mat.fgbg) <- N.mat.fgbg$gene.uniq; N.mat.fgbg$gene.uniq <- NULL
  
  # remove column sums == 0
  N.mat.fgbg[[names(which(colSums(N.mat.fgbg) == 0))]] <- NULL

  lda.out <- PenalizedLDA(N.mat.fgbg, labs, lambda = jlambda, K = K)
  return(lda.out)
}

BoxplotLdaOut <- function(out, jtitle = "Title"){
  out.df <- data.frame(xproj = out$xproj, lab = out$y)
  boxplot(xproj ~ lab, data = out.df, main = jtitle)
}

PlotLdaOut <- function(out, jtitle = "Title", jcex = 1){
  discrim.filt <- sort(out$discrim[which(out$discrim != 0)], decreasing = FALSE)
  cnames <- colnames(out$x)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
  textplot(x = seq(length(discrim.filt)), y = discrim.filt, words = cnames, main = jtitle, xlab = "Index", ylab = "Loadings", cex = jcex)
}

SortLda <- function(out){
  cnames <- colnames(out$x)[which(out$discrim != 0)][order(out$discrim[which(out$discrim != 0)], decreasing = FALSE)]
}

PlotSeparation <- function(out, jtitle){
  boxplot(list("Foreground" = out$xproj[out$y == 1], "Background" = out$xproj[out$y == 2]), main = jtitle, xlab = "Group 1: foreground. Group 2: background", ylab = "Projection")
}