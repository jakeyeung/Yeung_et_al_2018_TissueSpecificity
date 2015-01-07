LoadNormalizedArray <- function(normalized.array.path, remove.negs=TRUE){
  # Normalizing has caused strange behaviours in very flat genes, causing
  # the resulting output of one or two genes to be negative.
  # We need to remove these from subsequent analysis.
  
  normalized.array <- read.table(normalized.array.path)
  
  # How many have negative values? ------------------------------------------
  negs <- apply(normalized.array, 1, function(x){
    if (min(x) < 0){
      return(1)
    } else {
      return(0)
    }
  })
  
  problem.genes <- names(negs[which(negs == 1)])
  
  # Remove problem genes from analysis --------------------------------------
  genes <- rownames(normalized.array)
  filtered.genes <- genes[which(!genes %in% problem.genes)]
  normalized.array <- normalized.array[filtered.genes, ]
  
  print(paste("Removed problem genes:", problem.genes))
  return(normalized.array)
}

LoadRnaSeq <- function(rna.seq.path, handle.duplicates=TRUE){
  # tricky duplicate rownames. We'll fix that though.
  
  rna.seq.exprs <- read.table(rna.seq.path, header=TRUE, sep='\t')
  
  # Handle duplicate rownames: RNASEQ ------------------------------------
  
  rownames(rna.seq.exprs) <- make.names(rna.seq.exprs$gene, unique=TRUE)
  
  drop.cols <- c("gene")
  rna.seq.exprs <- rna.seq.exprs[, !(names(rna.seq.exprs) %in% drop.cols)]
  return(rna.seq.exprs)
}