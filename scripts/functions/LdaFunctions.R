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

PlotLdaOut <- function(lda.out, labs){
  plot(lda.out)
  text(x = seq(length(labs)), y = lda.out$discrim, labels = labs)
}