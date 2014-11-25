PlotRnaMicroarrayFit <- function(tissue, gene, coeff.mat, array.exprs, rna.seq.exprs, rna.tissue){
  # Plotting function when fitting array with expression RNA Seq with microarray.
  # 
  # Args:
  # tissue: string representing tissue name in ARRAY e.g. "Liver"
  # gene: gene to plot
  # coeff.mat: contains intercept and slope for each tissue for every gene. Generated from
  # fit_array_with_exprs.rna.seq.vs.array.R
  # 
  # array.exprs: expression of array
  # 
  # rna.seq.exprs: expression in rnaseq
  # 
  # rna.tissue: string representing tissue name in rna.seq .eg. "Liv"
  # 
  # 
  
  # plot for one tissue only
  t.grep <- tissue
  t.grep.rnaseq <- rna.tissue
  
  intercept <- coeff.mat[gene, paste0(tissue, '_intercept')]
  slope <- coeff.mat[gene, paste0(tissue, '_slope')]
  
  x <- as.matrix(array.exprs[gene, which(grepl(t.grep, colnames(array.exprs)))])
  y <- as.matrix(rna.seq.exprs[gene, which(grepl(t.grep.rnaseq, colnames(rna.seq.exprs)))])
  plot(x, y, xlab='Array exprs',
       ylab='RNA exprs',
       main=paste(gene, 'Intercept=', signif(intercept, digits=3), 'Slope=', signif(slope, digits=3)))
  abline(intercept, slope, lty='dotted')
}

PlotFitDiagnostics <- function(array.exprs.full,
                               array.exprs.subset,
                               rna.seq.exprs, 
                               clockgenes, 
                               coeff.mat,
                               outpath, 
                               tissue, 
                               tissue.rna.seq,
                               allow.negs=FALSE){
  if (missing(tissue.rna.seq) == TRUE){
    tissue.rna.seq <- tissue
  }
    
  pdf(outpath)
  par(mfrow=c(2, 2))
  for (gene in clockgenes){
    slope <- coeff.mat[gene, paste0(tissue, "_slope")]
    intercept <- coeff.mat[gene, paste0(tissue, "_intercept")]
    
    # get microarray and rnaseq
    x <- array.exprs.full[gene, grepl(tissue, colnames(array.exprs))]
    y.rna.seq <- rna.seq.exprs[gene, grepl(tissue.rna.seq, colnames(rna.seq.exprs))]
    # get indices of samples with both microarray and rnaseq
    x.rna.seq.i <- colnames(array.exprs.subset[gene, grepl(tissue, colnames(array.exprs.subset))])  # subset only ones that have mRNA matched
    # convert to indices
    x.rna.seq.i <- names(x) %in% x.rna.seq.i  # logical True/False
    # get subset x, containing microarray and rnaseq
    x.subset <- x[x.rna.seq.i]
    
    plot.symbols <- sapply(x.rna.seq.i, function(true.false){
      if (true.false == TRUE){
        symbol <- 1
      } else {
        symbol <- 8
      }
      })
    y <- slope * x + intercept
    if (allow.negs == FALSE){
      y[which(y < 1)] <- 1
    }
    plot(2^x, 2^y, main=paste("o=samp w/ RNASeq+array", gene), pch=c(plot.symbols), xlab="Observed microarray (normal)", ylab="Predicted expression (normal)")
    abline(h=0)
    
    plot(x, y, main=paste("log2", gene, tissue), pch=c(plot.symbols), xlab="Observed microarray (log2)", ylab="Predicted expression (log2)")
    abline(h=0)
    
    # Plot raw on normal scale
    
    plot(2^x.subset, 2^y.rna.seq, main=paste("RNASeq vs Array: normal scale"), xlab="Microarray normal scale", ylab="DESeq normalized counts")
    
    # Plot raw on log2 scale
    plot(x.subset, y.rna.seq, 
         main=paste(gene, 'Intercept=', signif(intercept, digits=3), 
                    'Slope=', signif(slope, digits=3)), 
         xlab="Microarray log2 scale", 
         ylab="DESeq normalized counts (log2)") 
    abline(intercept, slope, lty='dotted')
    }
  dev.off() 
}