# 2016-08-16
# Figures for paper: include hogenesch, liver kidney WTKO, nuclear proteomics, 4c-seq
# Jake Yeung
# 
# Other links to scripts not shown here:
# /home/yeung/projects/tissue-specificity/scripts/pca/pca_adjusted_microarray.label_variance.for_paper.R
# /home/yeung/projects/tissue-specificity/scripts/fourier/total_variance.noise_floor.hogenesch_and_liverWTKO.R
# /home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/downstream.three_way_contingency.R


rm(list=ls())
start <- Sys.time()

library(ggplot2)
library(PMA)
# detach("package:dplyr", unload=TRUE)  # sometimes necessary to solve some strange issues with PMA and dplyr
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
# source("scripts/functions/AnalyzeGeneEnrichment.R")
source("scripts/functions/FourierFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/ModelStrToModel.R")
source("scripts/functions/ProteomicsFunctions.R")
source("scripts/functions/CosSineFunctions.R")


# Inits -------------------------------------------------------------------

remove.kidney.outliers <- TRUE
remove.wfat <- TRUE
plot.i <- 1
tissue.order <- c('Liver','BFAT','Kidney','Lung','Adr','Mus','Heart','Aorta','Hypo','Cere','BS')

plot.dir <- "/home/yeung/projects/tissue-specificity/plots/primetime_plots_full_paper"
dir.create(plot.dir, showWarnings = FALSE)


tfs <- GetTFs(get.mat.only = TRUE)

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj", v=T)  # LivKid
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch
# model selection with g1000
load("Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj", v=T)  # use this, looks same as fits.best should be OK?
load("Robjs/dat.complex.fixed_rik_genes.Robj", v=T)
if (remove.wfat){
  dat.complex <- subset(dat.complex, tissue != "WFAT")
  dat.long <- subset(dat.long, tissue != "WFAT")
}

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

# load proteomics

prot.long <- LoadProteomicsData()

# 
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", v=T)

# Do tables for paper -----------------------------------------------------

maindir <- "/home/yeung/projects/tissue-specificity/tables"

# write microarray-RNASeq merged table

outf <- "1_Hogenesch-Microarray-RNASeq-Merged.txt"
write.table(dat.long, file = file.path(maindir, outf), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# write RNA-Seq WT-KO table

outf2 <- "2_Liver-Kidney-WTKO-RNASeq.txt"
write.table(subset(dat.wtko, !is.na(gene)), file = file.path(maindir, outf2), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# write fits for Hogenesch
outf3 <- "3_Hogenesch-Model-Selection.txt"
fits.best.txtout <- fits.best; fits.best.txtout$param.list <- NULL
fits.best.txtout$param.list <- sapply(fits.best$param.list, function(l){
  # Make a long string with Paramname = Value
  outstr <- paste(names(l), as.character(l), sep = "=", collapse = ";")
})



RenameCnameStr <- function(dat, old, new){
  return(dplyr::rename_(dat, .dots=setNames(list(old), new)))
}

fits.best.txtout <- RenameCnameStr(fits.best.txtout, "model", "Model Name (distinct rhythms are separated by semicolon, while shared rhythms separated by comma)")
fits.best.txtout <- RenameCnameStr(fits.best.txtout, "weight", "Posterior Probability")
fits.best.txtout <- RenameCnameStr(fits.best.txtout, "param.list", "Parameter List (Each name=value pair separated by semicolon. Tissue names denote intercept estimates for RNA-Seq or microarray. Phase-amp denotes rhythmic parameters for tissues. Shared tissue rhythms separated by comma.)")
fits.best.txtout$weight.raw <- NULL

write.table(fits.best.txtout, file = file.path(maindir, outf3), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# write fits for WT KO 
outf4 <- "4_LivKidWTKO-Model-Selection.txt"


fits.long.filt.txtout <- fits.long.filt; fits.long.filt.txtout$param.list <- NULL
fits.long.filt.txtout$param.list <- sapply(fits.long.filt$param.list, function(l){
  outstr <- paste(names(l), as.character(l), sep = "=", collapse = ";")
})
fits.long.filt.txtout <- RenameCnameStr(fits.long.filt.txtout, "model", "Model Name (distinct rhythms are separated by semicolon, while shared rhythms separated by comma)")
fits.long.filt.txtout <- RenameCnameStr(fits.long.filt.txtout, "weight", "Posterior Probability")
fits.long.filt.txtout <- RenameCnameStr(fits.long.filt.txtout, "param.list", "Parameter List (Each name=value pair separated by semicolon). Tissue names denote intercept estimates for RNA-Seq or microarray. Phase-amp denotes rhythmic parameters for tissues. Shared tissue rhythms separated by comma.)")
fits.long.filt.txtout$weight.raw <- NULL

write.table(fits.long.filt.txtout, file = file.path(maindir, outf4), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

