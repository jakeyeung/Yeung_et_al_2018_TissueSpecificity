# Functions involving linear regressions and setting up matrices to get it done.
# used in fit_array_with_exprs.rnaseq.vs.array.R and other similar scripts.
# November 24 2014
# Jake Yeung

InitCoeffMat <- function(row.names, col.names){
  # Prepare empty matrix of NAs with rownames as supplied.
  # Each colname is split into two columns in coeff.mat, a samp1_intercept and samp1_slope
  # 
  # Args
  # row.names: Rownames of coeff.mat
  # colnames: Colnames of coeff.mat
  # 
  # Output
  # coeff.mat: empty matrix of NAs containing rownames and colnames.
  
  coeff.mat <- matrix(NA, nrow=length(row.names), ncol=(length(col.names) * 2))
  # set up row and colnames of coeff.mat
  rownames(coeff.mat) <- row.names
  coeff.mat.colnames <- rep(NA, 2 * length(col.names))
  for (i in 1:length(col.names)){
    col <- col.names[i]
    coeff.mat.colnames[i * 2] <- paste0(col, '_slope')
    coeff.mat.colnames[i * 2 - 1] <- paste0(col, '_intercept')
  }
  colnames(coeff.mat) <- coeff.mat.colnames
  return(coeff.mat)
}

LmGeneTissue <- function(array.subset, rna.seq, row.names, 
                         tissue.names, n.samps){
  # Do lm fit to each tissue, for each gene.
  # 
  # Args:
  # array.subset: exprs of array, has same colnames as rna.seq
  # rna.seq: exprs of rnaseq
  # row.names: gene names
  # tissue.names: Tissue names
  # n.samps: samples in each tissues
  
  coeff.mat <- InitCoeffMat(row.names, tissue.names)
  
  for (gene in row.names){
    for (tissue.i in 1:length(tissue.names)){
      tissue <- tissue.names[tissue.i]
      tissue.i.end <- tissue.i * n.samps
      tissue.i.start <- tissue.i.end - n.samps + 1
      y <- rna.seq[gene, tissue.i.start:tissue.i.end]
      x <- array.subset[gene, tissue.i.start:tissue.i.end]
      fit <- lm(y ~ x)
      intercept <- fit$coefficient[1]
      slope <- fit$coefficient[2]
      coeff.mat[gene, paste0(tissue, '_slope')] <- slope
      coeff.mat[gene, paste0(tissue, '_intercept')] <- intercept
    }
  }
  return(coeff.mat)
}


AdjustArrayToRnaSeq <- function(array.exprs, coeff.mat, tissue.names){
  # init output df
  array.exprs.adjusted <- matrix(NA, nrow=nrow(array.exprs), ncol=ncol(array.exprs),
                                 dimnames=list(rownames(array.exprs),
                                               colnames(array.exprs)))
  for (tissue in tissue.names){
    intercept <- coeff.mat[, paste0(tissue, '_intercept')]
    slope <- coeff.mat[, paste0(tissue, '_slope')]
    tissue.exprs.array <- array.exprs[, grepl(tissue, 
                                              colnames(array.exprs))]
    tissue.exprs.array.normalized <- intercept + slope * tissue.exprs.array
    # write to adjusted exprs
    array.exprs.adjusted[, grepl(tissue, 
                                 colnames(array.exprs.adjusted))] <- tissue.exprs.array.normalized
  }
  return(array.exprs.adjusted)
}

