# 2016-08-05 
# Jake Yeung
# compare_hogenesch_tissue_wide_in_livkid_wtko.R

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/PlotFunctions.R")
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
source("scripts/functions/ProteomicsFunctions.R")

# Load hogenesch data -----------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
load("Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == "g=1001")
fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)

# Get tissue wide genes ---------------------------------------------------

genes.tw <- as.character(subset(fits.long, n.rhyth >= 8)$gene)

fits.hog <- subset(fits.long.filt, gene %in% genes.tw)


# Count models ------------------------------------------------------------

fits.sum <- fits.hog %>%
  group_by(model) %>%
  summarise(count = length(gene)) %>%
  arrange(desc(count))

fits.sum <- OrderDecreasing(fits.sum, jfactor = "model", jval = "count")
ggplot(fits.sum, aes(x = model, y = count)) + geom_bar(stat = "identity")
