Transcript2Gene <- function(gene.list, return.original=TRUE) {
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list,attributes=c("ensembl_transcript_id", "external_gene_name"),
               filters=c("ensembl_transcript_id"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}

AppendGeneID <- function(dat){
  genes.tens <- rownames(dat)
  Gene.ID <- Transcript2Gene(genes.tens)
  dat <- cbind(Gene.ID, dat)
  return(dat)
}
