# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")

library(topGO)
library(org.Mm.eg.db)


# Functions ---------------------------------------------------------------

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

# Create gene symbol to entrez object -------------------------------------

sym2entrez <- as.list(org.Mm.egALIAS2EG)
sym2entrez <- sym2entrez[!is.na(sym2entrez)]


# Create entrez to GO ID object -------------------------------------------

entrez2GO.obj <- org.Mm.egGO
mapped.genes <- mappedkeys(entrez2GO.obj)
entrez2GO.full <- as.list(entrez2GO.obj[mapped.genes])
# collapse sublist of GO-IDs to vector
entrez2GO <- lapply(entrez2GO.full, function(x){
  return(names(x))  # lapply is nice
})


# geneID2GO<- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))


# Run enrichment ----------------------------------------------------------

top.genes <- read.table("plots/nconds/7_conds_filtered/7_conds_filtered4.txt")
top.genes <- as.character(unlist(top.genes))  # to get correct naming after conversion

all.genes <- rownames(dat.with.fit.filtered)  # from run_nconds.R

genes.entrez <- ConvertSym2Entrez(all.genes, sym2entrez)
top.genes <- ConvertSym2Entrez(top.genes, sym2entrez)

sel.genes <- factor(as.integer(genes.entrez %in% top.genes))
names(sel.genes) <- genes.entrez

GOdata <- new("topGOdata",
              description = "Test",
              ontology = "BP",
              allGenes = sel.genes,
              nodeSize = 10,
              annot = annFUN.gene2GO,
              gene2GO = entrez2GO)