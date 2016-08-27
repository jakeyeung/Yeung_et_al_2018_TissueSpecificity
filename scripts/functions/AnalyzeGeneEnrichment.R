# AnalyzeGeneEnrichment.R
# using topGO
# 25 Feb 2014
#
# if necessary:
# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")


library(topGO)
library(org.Mm.eg.db)

CreateSym2Entrez <- function(){
  sym2entrez <- as.list(org.Mm.egALIAS2EG)
  sym2entrez <- sym2entrez[!is.na(sym2entrez)]
  return(sym2entrez)
}

CreateEntrez2GO <- function(){
  entrez2GO.obj <- org.Mm.egGO
  mapped.genes <- mappedkeys(entrez2GO.obj)
  entrez2GO.full <- as.list(entrez2GO.obj[mapped.genes])
  # collapse sublist of GO-IDs to vector
  entrez2GO <- lapply(entrez2GO.full, function(x){
    return(names(x))  # lapply is nice
  })
  return(entrez2GO)
}

ConvertSym2Entrez <- function(gene.list, sym2entrez){
  
  # convert to entrez, auto filtering here
  gene.list.entrez <- sapply(gene.list, function(g){
    entrez <- sym2entrez[[g]]
    if (!is.null(entrez)){
      return(entrez)
    }
  }, simplify = TRUE)
  gene.list.entrez <- gene.list.entrez[!sapply(gene.list.entrez, is.null)]
  return(unlist(gene.list.entrez))
}

AnalyzeGeneEnrichment <- function(genes.bg, genes.hit, 
                                  sym2entrez,
                                  entrez2GO,
                                  convert.sym.to.entrez = TRUE,
                                  which.ontology = "BP", 
                                  write.path = FALSE,
                                  node.size = 5, 
                                  FDR.cutoff = 0.05,
                                  return.GOdata = FALSE){
  # Analyze gene enrichment given background and hit genes.
  #
  # INPUT:
  # genes.bg: vector of background genes
  # genes.hit: vector of hit genes
  # sym2entrez: Provide a mapping from gene symbol to entrez, see example via CreateSym2Entrez(). If missing, creates it via CreateSym2Entrez()
  # entrez2GO: mapping from entrez to GO terms, example CreateEntrez2GO(). If missing, creates it via CreateEntrez2Go()
  # convert.sym.to.entrez: convert gene symbols to entrez ID. Set to FALSE if genes already in entrez ID
  # write.path: path to write topGO results object. If FALSE, does not write to file
  # which.ontology: "BP", "MF", "CC"
  # node.size: prune GO hierarchy from terms which have less than node.size genes
  # FDR.cutoff: Fisher's exact test FDR-adjusted pvals cutoff to be significant
  # write.path: if false jsut return object, if true also write to file given by write.path
  #
  # OUTPUT:
  # all.res: topGO results object
  #
  if (missing(sym2entrez)){
    sym2entrez <- CreateSym2Entrez()
  }
  if (missing(entrez2GO)){
    entrez2GO <- CreateEntrez2GO()
  }
  
  if (convert.sym.to.entrez){
    genes.bg <- ConvertSym2Entrez(genes.bg, sym2entrez)
    genes.hit <- ConvertSym2Entrez(genes.hit, sym2entrez)
  }
  
  # Select genes in bg that are in hit, binary matrix.
  # used in topGO new() function
  sel.genes <- factor(as.integer(genes.bg %in% genes.hit))
  names(sel.genes) <- genes.bg
  
  # Get topGO object
  GOdata <- new("topGOdata",
                ontology = which.ontology,
                allGenes = sel.genes,
                nodeSize = node.size,
                annot = annFUN.gene2GO,
                gene2GO = entrez2GO)
  # return(GOdata)
  
  # Run enrichment ----------------------------------------------------------
  
  result.fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  n.tests <- length(score(result.fisher))
  
  # Look at results ---------------------------------------------------------
  
  all.res <- GenTable(GOdata, 
                      classicFisher = result.fisher,
                      orderBy = "classicFisher",
                      ranksOf = "classicFisher",
                      topNodes = n.tests)
  # adjust pvalues
  all.res$FDRadj <- p.adjust(all.res$classicFisher, method = "BH")
  
  
  # Filter table by pvalue --------------------------------------------------
  
  all.res <- all.res[all.res$FDRadj <= FDR.cutoff, ]
  print(paste(nrow(all.res), "GO.terms found enriched."))
  if (nrow(all.res) <= 1){
    print("No enrichment found below cutoff")
    return(all.res)  # just skip it and return an empty DF
  }
  
  # Add all genes that were considered (can recapitulate contTable) ---------
  N.genes <- length(genes(GOdata))
  
  
  all.res$N.genes <- N.genes
  
  # Write to file -----------------------------------------------------------
  if (write.path != FALSE){
    write.table(all.res, file = write.path, 
                quote = FALSE, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = TRUE) 
  }
  if (!return.GOdata){
    return(all.res)
  } else {
    # https://support.bioconductor.org/p/65856/#66677
    # add gene names to all.res
    all.res$genes <- sapply(all.res$GO.ID, function(x){
      genes <- genesInTerm(GOdata, x) 
      genes[[1]][genes[[1]] %in% genes.hit]
    })
    # convert to gene symbol
    entrez2sym <- as.list(org.Mm.egSYMBOL)
    entrez2sym <- entrez2sym[!is.na(entrez2sym)]
    all.res$genes <- sapply(all.res$genes, function(g){
      return(unlist(entrez2sym[g], use.names = FALSE))
    })
    return(all.res)
    # return(list(res=all.res, godata=GOdata, genes=genes.hit))
  }
}




EnrichmentBinnedToFile <- function(sorted.hits, bin.vector, fname.base, sym2entrez, entrez2GO){
  # Take sorted hits and take top bin.vector[i] genes for enrichment
  # to see how enrichment of certain processes evolve over time
  # 
  # sorted.hits: list of genes sorted by "relevance" 
  # bin.vector: vector of integers to bin genes (take top 10 genes, then top 20 ... etc )
  # fname.base: full path of filename, we will add _i_Ontology.GOtop to end of file
  
  for (i in bin.vector){
    genes.hit <- sorted.hits[1:i]
    for (onto in c("BP", "MF", "CC")){
      fname.out <- paste0(fname.base, '_', i, '_', onto, '.GOtop')
      res <- AnalyzeGeneEnrichment(sorted.hits, genes.hit, 
                                   sym2entrez, entrez2GO, 
                                   convert.sym.to.entrez = TRUE, 
                                   which.ontology = onto, 
                                   write.path = fname.out, 
                                   node.size = 5, 
                                   FDR.cutoff = 0.05)
    }
  }
}

GetGoObject <- function(genes.hit, genes.bg, onto, node.size, sym2entrez, entrez2GO){
  # Run GO
  
}

GetEnrichment <- function(sorted.hits, bin, onto, sym2entrez, entrez2GO){
  # Same as getEnrichmentOverBins but for just one bin, that way we can parallelize it 
  genes.hit <- sorted.hits[1:bin]  
  res <- AnalyzeGeneEnrichment(sorted.hits, genes.hit, 
                               sym2entrez, entrez2GO, 
                               convert.sym.to.entrez = TRUE, 
                               which.ontology = onto, 
                               write.path = FALSE, 
                               node.size = 5, 
                               FDR.cutoff = 1)
  return(res)
}

GetEnrichmentParallel <- function(sorted.hits, bin.vector, onto, sym2entrez, entrez2GO, n.cores = 4){
  # Run GetEnrichment in parallel over all bins in bin.vector
  library(parallel)
  print("Running GO (~4 minutes)")
  start <- Sys.time()
  res.split <- mclapply(bin.vector, function(bin){
    genes.hit <- sorted.hits[1:bin]
    res <- AnalyzeGeneEnrichment(sorted.hits, genes.hit, 
                                 sym2entrez, entrez2GO, 
                                 convert.sym.to.entrez = TRUE, 
                                 which.ontology = onto, 
                                 write.path = FALSE, 
                                 node.size = 5, 
                                 FDR.cutoff = 1)
    res$bin <- bin
    return(res)
  }, mc.cores = n.cores)
  print(Sys.time() - start)
  res.all <- do.call(rbind, res.split)
  return(res.all)
}

GetEnrichmentOverBins <- function(sorted.hits, bin.vector, onto, go.terms, sym2entrez, entrez2GO){
  # Iterate over bins and find enrichment for a specific go term
  # track this enrichment as we sweep across the bins
  
  # init output matrix: 8 columns ("GO.ID", "Term", "Annotated", "Significant", "Expected", "classicFisher", "FDRadj", "N.genes", "bin")
  res.out <- matrix(NA, nrow = length(bin.vector) * length(go.terms), ncol = 9)
  colnames(res.out) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "classicFisher", "FDRadj", "N.genes", "bin")
  
  rowcount <- 1
  for (i in 1:length(bin.vector)){
    n.genes <- bin.vector[i]
    
    genes.hit <- sorted.hits[1:n.genes]
    
    res <- AnalyzeGeneEnrichment(sorted.hits, genes.hit, 
                                 sym2entrez, entrez2GO, 
                                 convert.sym.to.entrez = TRUE, 
                                 which.ontology = onto, 
                                 write.path = FALSE, 
                                 node.size = 5, 
                                 FDR.cutoff = 1)
    res.sub <- subset(res, Term %in% go.terms)
    res.sub$bin <- n.genes
    # make into matrix, makes life easier later
    res.sub <- as.matrix(res.sub)
    res.sub[, 3:ncol(res.sub)] <- as.numeric(res.sub[, 3:ncol(res.sub)])
    
    rowstart <- rowcount
    rowend <- rowstart + length(go.terms) - 1
    res.out[rowstart:rowend, ] <- res.sub
    rowcount <- rowcount + length(go.terms)
  }
#     # make to dataframe and unfactor things
#   res.out <- data.frame(noquote(res.out))
#   for (cname in c("Annotated", "Significant", "Expected", "classicFisher", "FDRadj")){
#     res.out[[cname]] <- as.numeric(res.out[[cname]])
#   }
#   res.out <- data.frame(res.out)
#   for (cname in c("Annotated", "Significant", "Expected", "classicFisher", "FDRadj", "bin")){
#     res.out[[cname]] <- as.numeric(as.character(res.out[[cname]]))
#   }
  return(res.out)
}