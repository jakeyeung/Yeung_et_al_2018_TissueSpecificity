LoadSitecounts <- function(N.path, N.promoterpath){
  if (missing(N.path)){
    N.path <- "data/sitecounts/motevo/sitecount_matrix_geneids"
  } 
  if (missing(N.promoterpath)){
    N.promoterpath <- "data/sitecounts/motevo/mara_promoters_gene_name_association.bed"
  }
  N <- read.table(N.path, header = TRUE, sep = '\t')
  rnames <- make.names(N$Gene.ID, unique = TRUE)
  rownames(N) <- rnames
  N$Gene.ID <- NULL
  
  N.promoter <- read.table(N.promoterpath, header = FALSE, sep = '\t')
  N.promoter <- data.frame(Gene.ID=rnames, maraid=N.promoter$V4, saeedid=N.promoter$V10)
  rownames(N.promoter) <- rnames
  return(list(N=N, N.promoter=N.promoter))
}