# 2016-08-30
# Jake YEung
# go_terms_analysis_with_phase_amplitude.R
# redo but load ALL terms over time, then not by "removing IDs" but by "choosing IDs"

rm(list=ls())
start <- Sys.time()

textsize <- 5
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
jtiss.lst <- list(list("a"=c("Liver_SV129,Liver_BmalKO"), 
                       # "b"=c("GO:0006270", "GO:0042254", "GO:0051917"),  # DNA rep, ribi, fibrinolysis
                       # "b"=c("GO:0006270", "GO:0042254", "GO:0032869"),  # DNA rep, ribi, response to insulin
                       "b"=c("GO:0006270", "GO:0042254", "GO:1900077"),  # DNA rep, ribi, neg response to insulin
                       # "b"=c("GO:1900077", strsplit("GO:0038014,GO:0043569,GO:0046627,GO:0046676,GO:0061179,GO:1900077,GO:1903105,GO:2001274", ",")[[1]]),  # neg response to insulin all 
                       # "b"=c("GO:0043434", "GO:0032868", "GO:0006111"), 
                       # "b"=c("GO:0051917", "GO:0006270"), 
                       "c"=4, 
                       "d"=c("Mcm2", "Mcm3", "Mcm4", "Mcm5", "Pml",   # DNA
                             "Egr1", "Kat2b", "Prkcd", "Pparg", "Cpeb2", "Ptprf", "Sorbs1",  # Response to insulin
                             "Erbb4", "Prkcd", "Fgf21",  # Positive reg
                             "Kat2b", "Sesn2", "Foxo1",  # Gluconeo
                             "Mafb", "Lpin1", "Jun", "Ddit3", "Cebpb", "Pik3r1", "Nop58", "Nop56", "Rpp38",  # ribi and misc
                             "Kank1", "Rps6kb1", "Grb10")),  # neg of insulin 
                  list("a"=c("Liver_SV129"),
                       # "b"=c("GO:0044262", "GO:0044255"),
                       # "b"=c("GO:0006793", "GO:0034637", "GO:0044255"),  # carbo, phospho, lipid,
                       "b"=c("GO:0034637", "GO:0009117", "GO:0044255"),  # replace phospho with puridine
                       # "b"=c("GO:0034637", "GO:0030301", "GO:0044255"),  # cholestrol
                       # "b"=c("GO:0071229","GO:0034637","GO:0006796","GO:0006793","GO:0000271","GO:0051865"),  # top guys
                        #"b"=c("GO:0006633", "GO:0006511", "GO:0055114", "GO:0070542", "GO:0007584", "GO:0009056", "GO:0006739", "GO:0006793", "GO:0051384"),
                       "c"=5,
                       "d"=c("Elovl3", "Insig2", "Adh4", "Hsd3b7", "Lrp5", "Sdr42e1", "Cyp8b1",
                             "Gck", "Pklr", "Ppp1r3b", "Colgalt2", 
                             "Upp2", "Rgs16", "Abcg5", "Por", "Lipg", "Tff3", "Fkbp4", "Nrg4",
                             "Slc45a3", "Pik3ap1", "Mreg", "Slc44a1", "Slc4a4")),
                  list("a"=c("Kidney_SV129"),
                       # "b"=c("GO:0071804", "GO:0090317"),
                       # "b"=c("GO:0098656", "GO:0006814"),
                       # "b"=c("GO:0098656","GO:0015727","GO:0035873","GO:0035879","GO:0015701"),
                       # "b"=c("GO:0051224","GO:0015711", "GO:0006814"),  # organic anion, neg reg of protein transport (DMNT not studied in kidney, ignore), sodium
                       "b"=c("GO:0015711", "GO:0006814"),  # organic anion,  sodium
                       # "b"=c("GO:0098656","GO:0034220","GO:0051224","GO:0055085","GO:0015711"),
                       # "b"=c("GO:0015727","GO:0035873","GO:0035879","GO:0015701"),
                       # "b"=c("GO:0006633", "GO:0006511", "GO:0055114", "GO:0070542", "GO:0007584", "GO:0009056", "GO:0000188", "GO:0048015", "GO:0006811", "GO:0055085"),
                       "c"=2.5,
                       "d"=c("Tfcp2", "Angpt1", "Igf1r",
                             "Slc12a6", "Prnp", "Lrrc52",
                             "Clcn2", "Spred3", "Dusp9",
                             "Slc16a1", "Slc7a8", "Slc39a5", "Clcn2", "Slc22a4", "Slc41a1", "Slc9a3", "Rhobtb1", "Trib2", "Slc6a4",
			                       "Slc7a8")),
                  list("a"=c("Liver_BmalKO"),
                       "b"=c("GO:0006695", "GO:0006695"),
                       "c"=5.5,
                       "d"=c("Akr1b7", "Gpam",
                             "Dhcr24", "Mvk",
                             "Ar", "Egfr", "Gstp1", "Traf6",
                             "Onecut1", "Uso1", "Rsg1", "Sec24d", "Gm9493", "Iars")))

jtiss.onto <- jtiss.lst[[3]]

plotdir <- "/home/yeung/projects/tissue-specificity/plots/liver_kidney_modules_with_GO"
pdf(file.path(plotdir, paste0("liv_kid_with_GO.phasewindow.rm_outliers.fixkidney.", remove.kidney.outliers, ".keepGOterms.neginsulin.pdf")))
  mclapply(jtiss.lst, function(jtiss.onto){
    # jtiss.onto <- jtiss.lst[[1]]  # c(jmodels, ontology), remove last element to get jmodels, last element is ontology
    source("scripts/functions/AnalyzeGeneEnrichment.R")
    comp <- 1
    plots <- expandingList()
    jmod.long <- jtiss.onto[1][[1]]
    GO.keep <- jtiss.onto[2][[1]]
    max.amp <- jtiss.onto[3][[1]]
    show.genes <- jtiss.onto[4][[1]]
    # Load GOenrichment results
    load(paste0("Robjs/GO_analysis/model", jmod.long, ".all.Robj"), v=T); enrichment <- CopyZT0ToZT24(enrichment, jorder = TRUE, convert.cname = TRUE)
    enrichment <- subset(enrichment, GO.ID %in% GO.keep)
    
    genes <- as.character(subset(fits.long.filt, model %in% jmod.long)$gene)
    dat.sub <- subset(dat.freq, gene %in% genes)
    s <- SvdOnComplex(dat.sub, value.var = "exprs.transformed")
    eigens <- GetEigens(s, period = 24, comp = comp, label.n = 0, eigenval = TRUE, adj.mag = TRUE, constant.amp = textsize, peak.to.trough = TRUE, label.gene = show.genes, dotsize = 3, dot.col = "gray85")
    # eigens <- GetEigens(s, period = 24, comp = comp, label.n = 50, eigenval = TRUE, adj.mag = TRUE, constant.amp = textsize, peak.to.trough = TRUE, label.gene = NA, dotsize = 3, dot.col = "gray85")
    
   
    eigentiss <- eigens$u.plot + ylab("ZT") + ggtitle("")
    eigengene <- eigens$v.plot + ylab("ZT") + xlab("Tissue Weights") + ggtitle("")
    
    cbPalette <- c("#E69F00", "#56B4E9", "#CC79A7", "#D55E00", "#009E73", "#0072B2", "#999999", "#F0E442")
    plots$add(eigentiss)
    plots$add(eigengene)
    plots$add(PlotGOTermsPhase(enrichment, "", cbPalette, to.append = eigentiss, phasestr = "phase", ampstr = "amp", amp.max = max.amp))
    plots$add(PlotGOTermsPhase(enrichment, "", cbPalette, to.append = eigentiss, phasestr = "phase", ampstr = "amp", amp.max = max.amp) + theme(legend.position = "right"))
    return(plots$as.list())
  }, mc.cores = length(jtiss.lst))
dev.off()
