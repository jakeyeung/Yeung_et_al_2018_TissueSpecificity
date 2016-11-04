# 2016-11-04
# Explore GO terms that plots trees

# source("http://www.bioconductor.org/biocLite.R")
# biocLite()
# biocLite("graph")
# biocLite("Rgraphviz") #before installing this package, the user has to install the graphviz 
# biocLite("SparseM")
# biocLite("GOFunction")

# Load --------------------------------------------------------------------

library("Rgraphviz")
library("GOFunction")
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(GO.db)
library(graph)

library(devtools)
# devtools::install_github("jakeyeung/GOFunction")
devtools::install_local("/home/yeung/projects/GOFunction")

source("scripts/functions/GraphFunctions.R")

data(exampledata)

sigTerm <- GOFunction(interestGenes, refGenes, organism="org.Hs.eg.db", ontology="BP", fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05, poth=0.05, peth=0.05, bmpSize=2000, filename="sigTerm")

organism="org.Hs.eg.db"; ontology="BP"; fdrmethod="BY"; fdrth=0.05; ppth=0.05; pcth=0.05; poth=0.05; peth=0.0; bmpSize=2000; filename="sigTerm"


# Hack code to get graph visualization without output to file -------------



## Extracting the gene annotation data
.sql <-  paste("select distinct t1.gene_id,t2.go_id",
               " from genes as t1 inner join",paste("go", tolower(ontology), "all", sep = "_"),
               " as t2 on t1._id=t2._id",seq="")
organism <- strsplit(organism,".db")
organism <- organism[[1]]
conn <- get(paste(organism, "_dbconn", sep = ""))()
generalAnn <- dbGetQuery(conn, .sql)
annRef <- generalAnn[generalAnn[,1] %in% refGenes,]
annInterest <- generalAnn[generalAnn[,1] %in% interestGenes,]

## Calculating the statistically significant terms for interestGenes
cat("Finding statistically significant terms...\n")
termInfo <- enrichmentFunction(annRef, annInterest,fdrmethod,fdrth)
sigTerm <- termInfo$sigTerm
if(nrow(sigTerm)==0){
  warning("There is no significant term! \n")
  return(NULL)
}
allTerm <- termInfo$allTerm

## Loading the GO structure data
conn <- get("GO_dbconn")()

## Finding the GO term name for statistically siginificant terms
.sql <- paste("select distinct go_id goid,term name from go_term where ontology='",
              toupper(ontology), "'", sep="")
allTermName <-  dbGetQuery(conn,.sql)
sigTermName <- allTermName[allTermName[,1] %in% sigTerm[,1],]
sigTermName <- sigTermName[order(sigTermName[,1]),]
sigTerm <- sigTerm[order(sigTerm[,1]),]
sigTerm$name <- sigTermName[,2]
sigTerm <- sigTerm[,c(1,6,2,3,4,5)]

## Finding the relationship between statistically significant terms
.sql <- paste("select distinct t1.go_id parentid,t2.go_id childid from ", paste("go", tolower(ontology),   "offspring", sep = "_"), " as t3 inner join  go_term as t1 on t1._id=t3._id inner join go_term as t2", " on t2._id=t3._offspring_id", sep="")
allTermRelation <-  dbGetQuery(conn,.sql)
sigTermRelation <- allTermRelation[(allTermRelation[,1] %in% sigTerm[,1]) & (allTermRelation[,2] %in% sigTerm[,1]),]
rm(allTermRelation)

## Finding the offspring terms for each statistically significant terms

## Local redundance
cat("Treating for local redundant terms...\n")
sigTerm_LocalRedun <- localRedundancy(sigTerm, generalAnn, sigTermRelation, annRef, annInterest, ppth, pcth)

##Global redundance

cat("Treating for global redundant terms...\n")
sigTerm_GlobalRedun <- globalRedundancy(generalAnn, sigTermRelation, annRef, annInterest, sigTerm_LocalRedun, poth,  peth)


## Display the GO DAG plot for the significant terms

cat("Visualizing the GO DAG...\n")

# take top 10 hits 
top.go.terms <- as.character(subset(sigTerm, pvalue < 1e-30)$goid)[1:2]
sigDAG <- createGODAG(top.go.terms, ontology)
allDAGTerm <- allTerm[allTerm[,1] %in% nodes(sigDAG),]

dagTermName <- allTermName[allTermName[,1] %in% allDAGTerm[,1],]
dagTermName <- dagTermName[order(dagTermName[,1]),]
allDAGTerm <- allDAGTerm[order(allDAGTerm[,1]),]
allDAGTerm$name <- dagTermName[,2]
allDAGTerm <- allDAGTerm[,c(1,6,2,3,4,5)]
sigTermID <- as.character(sigTerm[,1])
sigTerm_LocalRedunID <- as.character(sigTerm_LocalRedun[,1])
sigTerm_GlobalRedunID <- as.character(sigTerm_GlobalRedun[,1])

# BEGIN PLOTTING

graphAttrs <- getDefaultAttrs(layoutType = 'dot')
graphAttrs$cluster <- NULL
graphAttrs$node$shape <- 'ellipse'
graphAttrs$node$fontsize <- '20'

nodeAttrs <- list()
edgeAttrs <- list()

allTerm <- as.character(allDAGTerm[,1])
GO.names <- as.character(allDAGTerm$name)
nodeAttrs$label[allTerm] <- allTerm

rmLocalTerm <- setdiff(sigTermID, sigTerm_LocalRedunID)
nodeAttrs$color[rmLocalTerm] <- rep('red', length(rmLocalTerm))
nodeAttrs$shape[rmLocalTerm] <- rep('circle', length(rmLocalTerm))

rmGlobalTerm <- setdiff(sigTerm_LocalRedunID, sigTerm_GlobalRedunID)
nodeAttrs$color[rmGlobalTerm] <- rep('red', length(rmGlobalTerm))
nodeAttrs$shape[rmGlobalTerm] <- rep('box', length(rmGlobalTerm))
nodeAttrs$height[rmGlobalTerm] <- rep('0.7', length(rmGlobalTerm))
nodeAttrs$width[rmGlobalTerm] <- rep('0.7', length(rmGlobalTerm))

nodeAttrs$color[sigTerm_GlobalRedunID] <- rep('red', length(sigTerm_GlobalRedunID))
nodeAttrs$shape[sigTerm_GlobalRedunID] <- rep('rectangle', length(sigTerm_GlobalRedunID))
nodeAttrs$height[sigTerm_GlobalRedunID] <- rep('0.7', length(sigTerm_GlobalRedunID))
nodeAttrs$width[sigTerm_GlobalRedunID] <- rep('1.1', length(sigTerm_GlobalRedunID))

allDAGTerm[allDAGTerm[,5]<2.2E-16,5] <- 2.2E-16;
allDAGTerm[allDAGTerm[,6]<2.2E-16,6] <- 2.2E-16;
allDAGTerm$colorran <- round(log10(allDAGTerm[,6])-range(log10(allDAGTerm[,6]))[1] + 1)
mm <- max(allDAGTerm$colorran)
colorMap <- heat.colors(mm)
nodeAttrs$fillcolor[allTerm] <- unlist(lapply(allDAGTerm$colorran, function(x) return(colorMap[x])))

weightsList <- edgeWeights(sigDAG)
to <- lapply(weightsList, names)
from <- nodes(sigDAG)
edge.names <- paste(rep(from, listLen(to)), unlist(to), sep = "~")
edge.weights <- unlist(weightsList)
names(edge.weights) <- edge.names
##    0 for a is_a relation,  1 for a part_of relation
edgeAttrs$color <- ifelse(edge.weights == 0, 'black', 'red')

# nodeAttrs$label <- seq(length(nodeAttrs$label))
nodeAttrs$label[allTerm] <- GO.names
# nodeAttrs$cex[allTerm] <- 50
plot(sigDAG, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, main = "test")



