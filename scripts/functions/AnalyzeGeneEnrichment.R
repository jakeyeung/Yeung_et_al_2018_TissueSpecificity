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
                                  FDR.cutoff = 0.05){
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
  # write.path:
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
  names(sel.genes) <- genes.entrez
  
  # Get topGO object
  GOdata <- new("topGOdata",
                ontology = which.ontology,
                allGenes = sel.genes,
                nodeSize = node.size,
                annot = annFUN.gene2GO,
                gene2GO = entrez2GO)
  
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
  
  # Write to file -----------------------------------------------------------
  if (write.path != FALSE){
    write.table(all.res, file = write.path, 
                quote = FALSE, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = TRUE) 
  }
  return(all.res)
}