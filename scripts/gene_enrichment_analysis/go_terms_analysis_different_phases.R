# 2016-10-26
# Jake YEung
# go_terms_analysis_with_phase_amplitude.R
# copied from liver kidney WT KO analysis
# Create a phase-specific annotation (sliding window?) of functional enrichment. 

rm(list=ls())
start <- Sys.time()

textsize <- 7
dotsize <- 2
remove.kidney.outliers <- TRUE

library(ggplot2)
# detach("package:PMA", unload=TRUE)
# detach("package:plyr", unload=TRUE)
# detach("package:ggplot2", unload=TRUE)
# detach("package:reshape2", unload=TRUE)
library(dplyr)
library(parallel)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/FourierFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/ModelStrToModel.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/BiomartFunctions.R")


WriteListToFile <- function(lst, outf){
  sink(file = outf)
  for (g in lst){
    cat(g)
    cat("\n")
  }
  sink()
  return(NA)
}


Vectorize(IsBtwnTimes <- function(phase, tstart, tend){
  # check if phase (between 0 to 24, is between tstart and tend, considering the modulo)
  if (tend > tstart){
    # easy case
    is.btwn <- phase >= tstart & phase <= tend
  } else {
    # harder case, must consider the modulo
    is.btwn <- phase >= tstart | phase <= tend
  }
  # replace NAs with FALSE
  is.btwn[which(is.na(is.btwn))] <- FALSE
  return(is.btwn)
}, vectorize.args="phase")

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch

# Annotate each model with its major GO term ------------------------------

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)

fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)

# Remove kidney outliers (optional)
if (remove.kidney.outliers){
  # Kidney_SV129 genes contain some weird outliers, remove them
  outliers <- c("Cox8b", "Ucp1", "Cidea", "Flg", "Nr4a2")
  fits.long.filt <- subset(fits.long.filt, !gene %in% outliers)
  dat.wtko <- subset(dat.wtko, !gene %in% outliers)
}
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
dat.wtko.collapsed <- CollapseTissueGeno(dat.wtko)

# Plot for all genes in module --------------------------------------------

# from tissue_specificity_paper3.R
jmod1 <- "Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Liver_BmalKO.Kidney_SV129,Kidney_BmalKO"
jmod2 <- "Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129"
jtiss.lst <- list(c(ModelStrToModel(jmod1), "MF"),
                  c(ModelStrToModel(jmod2), "BP"),
                  c("Liver_SV129", "BP"), 
                  c("Liver_SV129,Liver_BmalKO", "BP"), 
                  c("Kidney_SV129", "MF"), 
                  c("Kidney_SV129,Kidney_BmalKO", "MF"), 
                  c("Liver_BmalKO", "BP"))

plots <- expandingList()
jmod.long <- "Liver_SV129,Liver_BmalKO"
jonto <- "BP"
comp <- 1
jcutoff <- 0.5
jtiss.onto <- paste(jmod.long, jonto)

genes <- as.character(subset(fits.long.filt, model %in% jmod.long)$gene)
dat.sub <- subset(dat.freq, gene %in% genes)
s <- SvdOnComplex(dat.sub, value.var = "exprs.transformed")
eigens <- GetEigens(s, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = textsize, peak.to.trough = TRUE, label.gene = c("Egr1", "Mafb", "Tfcp2"))

plots$add(eigens$u.plot + ylab("ZT") + ggtitle(""))
plots$add(eigens$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle(""))

genes.bg <- as.character(fits.long.filt$gene)

# write fg and bg to separate file to load into David:
genedir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/GO_analysis"

# GO terms for Liver WTKO module
GOterms <- c("GO:0006260", "GO:0042254", "GO:0032868", "GO:0043434", "GO:0046326")
# plot in 6 hour intervals
tstarts <- seq(0, 23)
enrichment <- lapply(tstarts, function(tstart){
  tend <- tstart + 6
  if (tend > 24){
    tend <- tend - 24
  }
  print(paste(tstart, tend))
  genes <- as.character(subset(fits.long.filt, model %in% jmod.long & IsBtwnTimes(phase.avg, tstart, tend))$gene)
  print(paste("Ngenes", length(genes)))
  genes.ensembl <- Gene2Ensembl(genes)
  enrichment <- GetGOEnrichment(genes.bg, genes, fdr.cutoff = 1, ontology = jonto, show.top.n = Inf, filter.GO.terms = GOterms)
  return(enrichment, !is.na(GO.ID))
})
