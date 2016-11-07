# 2016-11-06
# Jake Yeung
# load rhythmic GO terms and plot GO tree

rm(list=ls())

library(dplyr)
library(ggplot2)
library("Rgraphviz")
library("GOFunction")
library(org.Mm.eg.db)
library(GO.db)
library(graph)
library(ggplot2)
library(dplyr)
library(hash)

source("scripts/functions/GraphFunctions.R")
source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/HashFunctions.R")
source("scripts/functions/PhaseColorFunctions.R")

# Load rhythmic GO terms --------------------------------------------------

jmod <- "Kidney_SV129"
indir <- "/home/yeung/projects/tissue-specificity/Robjs/GO_analysis"
fname <- paste0("model", jmod, ".all.Robj")
inf <- file.path(file.path(indir, fname))

load(inf, v=T)

enrichment.sum <- enrichment %>%
  mutate(experiment = "GO", tissue = jmod, exprs = minuslogpval, time = tstart + 3) %>%   # define cols needed for FitRhythmic
  group_by(GO.ID, Term) %>%
  do(FitRhythmic(.)) %>%
  arrange(desc(amp))

# show top hits 
top.n <- 50
print(subset(data.frame(head(enrichment.sum, n = top.n)), select = c(GO.ID, Term, amp, phase, pval)))

# plot top hits
jGOID <- "GO:0090315"
jGOID <- "GO:0000001"
jtest <- subset(enrichment, GO.ID == jGOID)
jname <- jtest$Term
ggplot(jtest, aes(x = tstart + 3, y = minuslogpval)) + geom_point() + geom_line() + ggtitle(jname)


# Get DAG graph -----------------------------------------------------------

top.terms <- head(subset(data.frame(enrichment.sum %>% arrange(desc(amp))), select = c(GO.ID, Term, amp, phase, pval)), n = top.n)$GO.ID
top.terms <- c("GO:0006812", "GO:0006865")

ontology <- "BP"
sigDAG <- createGODAG(top.terms, ontology)

# Plot DAG graph ----------------------------------------------------------

pdf(paste0("plots/GO_analysis/GO_tree.", jmod, ".pdf"))
PlotDAGGraph(sigDAG, enrichment.sum, jfontsize = 75, jtitle = jmod, amp.cutoff = 0.3)
dev.off()
