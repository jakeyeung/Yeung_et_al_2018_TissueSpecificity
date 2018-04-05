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
source("scripts/functions/PlotFunctions.R")
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


load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch

# Annotate each model with its major GO term ------------------------------

jmeth <- "g=1001"
# load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
# load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
# load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)

fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)

# Remove kidney outliers (optional)
if (remove.kidney.outliers){
  # Kidney_SV129 genes contain some weird outliers, remove them
  outliers <- c("Cox8b", "Ucp1", "Cidea", "Flg", "Nr4a2")
  fits.long.filt <- subset(fits.long.filt, !gene %in% outliers)
  # dat.wtko <- subset(dat.wtko, !gene %in% outliers)
}
# dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
# dat.wtko.collapsed <- CollapseTissueGeno(dat.wtko)

# Plot for all genes in module --------------------------------------------

# from tissue_specificity_paper3.R
jmod.long <- c("Liver_SV129,Liver_BmalKO", "Liver_SV129")
jmod.str <- paste(jmod.long, collapse = "-")
jonto <- "BP"
comp <- 1
jcutoff <- 0.5
jtiss.onto <- paste(jmod.long, jonto)

genes.bg <- as.character(fits.long.filt$gene)

# write fg and bg to separate file to load into David:
genedir <- "/home/yeung/projects/tissue-specificity/data/gene_lists/GO_analysis"

# GO terms for Liver WTKO module
start <- Sys.time()
# DNA rep, response to insulin, ribosome biogenesis, glucose import ... 
# Liver WTKO GO terms
GOterms <- c("GO:0006260", "GO:0042254", "GO:0032868", "GO:0043434", "GO:0046326", "GO:0090526", "GO:0006111")  
# Liver WT GO terms 
GOterms <- c(GOterms, "GO:008152", "GO:0006633", "GO:0055114", "GO:46856", "GO:0070542", "GO:0006629", "GO:0006511", "GO:0051384", "GO:0007584", "GO:0005975", "GO:0006739", "GO:0009056")
# plot in 6 hour intervals
tstarts <- seq(0, 23)
enrichment <- mclapply(tstarts, function(tstart){
  source("scripts/functions/AnalyzeGeneEnrichment.R")
  tend <- tstart + 6
  if (tend > 24){
    tend <- tend - 24
  }
  print(paste(tstart, tend))
  genes <- as.character(subset(fits.long.filt, model %in% jmod.long & IsBtwnTimes(phase.avg, tstart, tend))$gene)
  print(paste("Ngenes", length(genes)))
  # genes.ensembl <- Gene2Ensembl(genes)
  enrichment <- GetGOEnrichment(genes.bg, genes, fdr.cutoff = 1, ontology = jonto, show.top.n = Inf, filter.GO.terms = GOterms)
  enrichment <- subset(enrichment, !is.na(GO.ID))
  enrichment$tstart <- tstart
  return(as.data.frame(subset(enrichment, !is.na(GO.ID))))
}, mc.cores = 12)
print(enrichment)
enrichment <- bind_rows(enrichment)
print(Sys.time() - start)

# save to output
fout <- paste0("/home/yeung/projects/tissue-specificity/Robjs/GO_analysis/model", jmod.str, ".Robj")
save(enrichment, file = fout)

# ggplot(subset(enrichment, GO.ID == "GO:0032868"), aes(x = tstart, y = minuslogpval)) + geom_line()
