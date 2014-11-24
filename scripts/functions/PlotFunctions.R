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