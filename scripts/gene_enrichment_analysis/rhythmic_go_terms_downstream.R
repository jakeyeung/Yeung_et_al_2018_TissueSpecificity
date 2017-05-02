# 2016-11-06
# Jake Yeung
# load rhythmic GO terms and plot GO tree

rm(list=ls())

setwd("~/projects/tissue-specificity")

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

manual.terms <- TRUE
manual.terms <- FALSE 

# Load rhythmic GO terms --------------------------------------------------

jmod <- "Liver_BmalKO"
jmod <- "Liver_SV129"
jmod <- "Kidney_SV129"
indir <- "/home/yeung/projects/tissue-specificity/Robjs/GO_analysis"
jmod <- "Liver_SV129,Liver_BmalKO"
jmod <- "Kidney_SV129,Kidney_BmalKO"

fname <- paste0("model", jmod, ".all.Robj")
inf <- file.path(file.path(indir, fname))

load(inf, v=T)

enrichment.sum <- enrichment %>%
  mutate(experiment = "GO", tissue = jmod, exprs = minuslogpval, time = tstart + 3) %>%   # define cols needed for FitRhythmic
  group_by(GO.ID, Term) %>%
  do(FitRhythmic(.)) %>%
  arrange(desc(amp))

# show top hits 
top.n <- 30
print(subset(data.frame(head(enrichment.sum, n = top.n)), select = c(GO.ID, Term, amp, phase, pval)))

# # plot top hits
# jGOID <- "GO:0090315"
# jGOID <- "GO:0000001"
jGOID <- "GO:0008643"
jGOID <- "GO:0044255"
jGOID <- "GO:0044262"
jGOID <- "GO:0015695"
jGOID <- "GO:0072337"
# search for a term
enrichment[grepl("cellular copper ion homeostasis", enrichment$Term), ]
enrichment[grepl("copper ion homeostasis", enrichment$Term), ]
enrichment[grepl("potassium ion transmembrane transport", enrichment$Term), ]
enrichment[grepl("cation transport", enrichment$Term), ]
jGOID <- "GO:0042254"  # ribi
jGOID <- "GO:0006487"  # N-glycosylation
jGOID <- "GO:0046326"  # glucose import
jGOID <- "GO:0006813"  # potassium ion transport
jGOID <- "GO:0043266"  # reg of pot ion transport
jGOID <- "GO:0044070"  # anion transport
jGOID <- "GO:0006878"  # cellular copper homeostasis
jGOID <- "GO:0055070"  # copper ion homeostasis 
jtest <- subset(enrichment, GO.ID == jGOID)
jname <- jtest$Term
ggplot(jtest, aes(x = tstart + 3, y = minuslogpval)) + geom_point() + geom_line() + ggtitle(jname) + ylab("Minus Log10 Pval") + xlab("ZT") + theme_bw()



# Get DAG graph -----------------------------------------------------------

if (!manual.terms){
  top.terms <- head(subset(data.frame(enrichment.sum %>% arrange(desc(amp))), select = c(GO.ID, Term, amp, phase, pval)), n = top.n)$GO.ID
  jfontsize <- 150
  jamp.cutoff <- 0.7
} else {
  if (jmod == "Liver_SV129"){
    jterms <- enrichment.sum[grepl("cellular carbohydrate metabolic|cellular lipid metabolic", enrichment.sum$Term), ]$GO.ID
    # top.terms <- c("GO:0044255", "GO:0044262")
    top.terms <- jterms
    jfontsize <- 50
    jamp.cutoff <- 0.5
  } else if (jmod == "Kidney_SV129"){
    jterms <- enrichment.sum[grepl("cation transport|amino acid transport", enrichment.sum$Term), ]$GO.ID
    jterms <- enrichment.sum[grepl("cation transport|amino acid transport|potassium ion transmembrane transport|cellular copper ion homeostasis", enrichment.sum$Term), ]$GO.ID
    top.terms <- jterms
    jfontsize <- 50
    jamp.cutoff <- 0.3
  } else if (jmod == "Liver_SV129,Liver_BmalKO"){
    jterms <- enrichment.sum[grepl("ribosome biogenesis|DNA replication initiation|response to insulin|regulation of fibrinolysis|positive regulation of glucose import", enrichment.sum$Term), ]$GO.ID
    # jterms <- enrichment.sum[grepl("ribosome biogenesis|DNA replication initiation|response to insulin|regulation of fibrinolysis|regulation of fibrinolysis", enrichment.sum$Term), ]$GO.ID
    top.terms <- jterms
    jamp.cutoff <- 0.7
    jfontsize <- 100
  } else if (jmod == "Liver_BmalKO"){
    top.terms <- c()
  }
}

ontology <- "BP"
sigDAG <- createGODAG(top.terms, ontology)

# Plot DAG graph ----------------------------------------------------------

pdf(paste0("plots/GO_analysis/GO_tree.", jmod, ".manualterms2.", manual.terms, ".pdf"))
PlotDAGGraph(sigDAG, enrichment.sum, jfontsize = jfontsize, jtitle = jmod, amp.cutoff = jamp.cutoff)
dev.off()
