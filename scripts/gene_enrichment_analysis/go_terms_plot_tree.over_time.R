# 2016-11-04
# Explore GO terms that plots trees

library("Rgraphviz")
library("GOFunction")
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GO.db)
library(graph)
library(ggplot2)
library(dplyr)
library(hash)

library(devtools)
# devtools::install_github("jakeyeung/GOFunction")
devtools::install_local("/home/yeung/projects/GOFunction")

source("scripts/functions/GraphFunctions.R")
source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/HashFunctions.R")
source("scripts/functions/PhaseColorFunctions.R")


# Load GO output ----------------------------------------------------------

inf <- "/home/yeung/projects/tissue-specificity/Robjs/GO_analysis/modelLiver_SV129,Liver_BmalKO.all.Robj"
load(inf, v=T)


# FInd oscillating  -------------------------------------------------------

enrichment.sum <- enrichment %>%
  mutate(experiment = "GO", tissue = "Liver", exprs = minuslogpval, time = tstart + 3) %>%   # define cols needed for FitRhythmic
  group_by(GO.ID, Term) %>%
  do(FitRhythmic(.))

# show top 50
head(subset(data.frame(enrichment.sum %>% arrange(desc(amp))), select = c(GO.ID, Term, amp, phase, pval)), n = 50)

# jtest <- subset(enrichment, Term == "ribosome biogenesis")
jGOID <- "GO:0009057"
jGOID <- "GO:0043434"
jGOID <- "GO:0031328"
jGOID <- "GO:0044260"
jGOID <- "GO:0006270"
jGOID <- "GO:0032868"
jtest <- subset(enrichment, GO.ID == jGOID)
jname <- jtest$Term
ggplot(jtest, aes(x = tstart + 3, y = minuslogpval)) + geom_point() + geom_line() + ggtitle(jname)

enrichment.sum$jlab <- apply(enrichment.sum, 1, function(row) ifelse(row[3] > 1, row[2], NA))

print(PlotAmpPhase(subset(enrichment.sum, amp > 1.25), lab_str = "jlab", textsize = 5) + xlim(c(0, 2)))


# Visualize top 25 GO terms in tree ---------------------------------------

top.terms <- head(subset(data.frame(enrichment.sum %>% arrange(desc(amp))), select = c(GO.ID, Term, amp, phase, pval)), n = 25)$GO.ID

ontology <- "BP"

## Finding the GO term name for statistically siginificant terms
conn <- get("GO_dbconn")()
.sql <- paste("select distinct go_id goid,term name from go_term where ontology='",
              toupper(ontology), "'", sep="")
allTermName <-  dbGetQuery(conn,.sql)
sigTermName <- allTermName[allTermName[,1] %in% top.terms,]
sigTermName <- sigTermName[order(sigTermName[,1]),]
sigTerm <- list()
sigTerm$name <- sigTermName[,2]

## Finding the relationship between statistically significant terms
.sql <- paste("select distinct t1.go_id parentid,t2.go_id childid from ", paste("go", tolower(ontology),   "offspring", sep = "_"), " as t3 inner join  go_term as t1 on t1._id=t3._id inner join go_term as t2", " on t2._id=t3._offspring_id", sep="")
allTermRelation <-  dbGetQuery(conn,.sql)
sigTermRelation <- allTermRelation[(allTermRelation[,1] %in% top.terms) & (allTermRelation[,2] %in% top.terms),]
rm(allTermRelation)

## Finding the offspring terms for each statistically significant terms

## Local redundance
# cat("Treating for local redundant terms...\n")
# sigTerm_LocalRedun <- localRedundancy(sigTerm, generalAnn, sigTermRelation, annRef, annInterest, ppth, pcth)

##Global redundance
# cat("Treating for global redundant terms...\n")
# sigTerm_GlobalRedun <- globalRedundancy(generalAnn, sigTermRelation, annRef, annInterest, sigTerm_LocalRedun, poth,  peth)


## Display the GO DAG plot for the significant terms

cat("Visualizing the GO DAG...\n")

# take top 10 hits 
pvalcutoff <- 1
sigDAG <- createGODAG(top.terms, ontology)



allDAGTerm <- allTermName$goid[allTermName$goid %in% nodes(sigDAG)]

dagTermName <- allTermName$goid[allTermName$goid %in% allDAGTerm]
dagTermName <- dagTermName[order(dagTermName)]
allDAGTerm <- allDAGTerm[order(allDAGTerm)]
allDAGTerm.name <- allTermName$name[allTermName$goid %in% allDAGTerm]

# allDAGTerm <- allDAGTerm[,c(1,6,2,3,4,5)]
sigTermID <- as.character(top.terms)

# BEGIN PLOTTING

graphAttrs <- getDefaultAttrs(layoutType = 'dot')
graphAttrs$cluster <- NULL
graphAttrs$node$shape <- 'ellipse'
graphAttrs$node$fontsize <- '20'

nodeAttrs <- list()
edgeAttrs <- list()

allTerm <- as.character(allDAGTerm)
GO.names <- as.character(allDAGTerm.name)
nodeAttrs$label <- nodes(sigDAG)

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
nodeAttrs$fixedsize[allTerm] <- TRUE
nodeAttrs$fontsize[allTerm] <- 300
# nodeAttrs$cex[allTerm] <- 50
# add node colors based on HSV
allTerm.df <- data.frame(goid = allDAGTerm, name = allDAGTerm.name, stringsAsFactors = FALSE)
pval.hash <- hash(enrichment.sum$GO.ID, enrichment.sum$pval)
amp.hash <- hash(enrichment.sum$GO.ID, enrichment.sum$amp)
phase.hash <- hash(enrichment.sum$GO.ID, enrichment.sum$phase)
allTerm.df$pval <- sapply(allTerm.df$goid, function(g) AssignHash(g, pval.hash, null.fill = NA))
allTerm.df$amp <- sapply(allTerm.df$goid, function(g) AssignHash(g, amp.hash, null.fill = NA))
allTerm.df$phase <- sapply(allTerm.df$goid, function(g) AssignHash(g, phase.hash, null.fill = NA))
allTerm.df$Color <- PhaseAmpPvalToColor(phase = allTerm.df$phase, amp = allTerm.df$amp, pval = allTerm.df$pval, 
                                        rotate.hr = -8, amp.k = 1, pval.k = Inf, method = "cutoff", black.to.white = TRUE)
color.hash <- hash(allTerm.df$goid, allTerm.df$Color)
nodeAttrs$fillcolor <- sapply(names(nodeAttrs$label), function(l) AssignHash(l, color.hash, null.fill = NA))
plot(sigDAG, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, main = "Testing")
