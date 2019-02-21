# Jake Yeung
# Date of Creation: 2019-02-21
# File: ~/projects/tissue-specificity/scripts/primetime_figures/create_primetime_objects.R
# Create primetime objects


rm(list=ls())
start <- Sys.time()

library(ggplot2)
# library(PMA)
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
# source("scripts/functions/LdaFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/GetTopMotifs.R")
source("scripts/functions/NcondsAnalysisFunctions.R")
source("scripts/functions/ModelStrToModel.R")
source("scripts/functions/ProteomicsFunctions.R")
source("scripts/functions/CosSineFunctions.R")
source("scripts/functions/HalfLifeFunctions.R")

dotsize <- 6
zscore.min <- 1.25
omega <- 2 * pi / 24
n <- 4  # amplitude scaling 
# mrna.hl <- 2.65  # shift TF activity by half-life

hr.shift <- 3
mrna.hl <- log(2) / (omega / tan(omega * hr.shift))  # convert hour shift to half-life 

# Functions ---------------------------------------------------------------

GetAmpPhaseFromActivities <- function(act.l, mrna.hl, jtiss = "Liver", jgeno = "WT"){
  act.l.complex <- ProjectWithZscore(act.l, omega, n = 4)
  act.l.complex$phase <- AdjustPhase(act.l.complex$phase, half.life = mrna.hl, fw.bw = "bw")
  return(act.l.complex)
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  # http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

RunPldaSystemsClock <- function(motifs, genes, clock.sys.gene.hash, N.long, jlambda = 0.035, take.max=FALSE){
  N.sub <- subset(N.long, motif2 %in% motifs & gene %in% genes)
  if (take.max){
    N.sub <- N.sub %>%
      group_by(motif, motif2, gene) %>%
      summarise(sitecount = max(sitecount))  # don't double count
    if (nrow(N.sub) == 1){
      warning("dplyr output has one row. Reload dplyr package")
    }
  }
  N.sub$clksys <- sapply(as.character(N.sub$gene), function(g) clock.sys.gene.hash[[g]])
  
  M.full <- dcast(N.sub, formula = gene + clksys ~ motif2, value.var = "sitecount", fun.aggregate = sum)
  
  M <- subset(M.full, select = c(-gene, -clksys)); rownames(M) <- M.full$gene
  M.labs <- as.numeric(as.factor(M.full$clksys))
  
  # END LOAD
  # load("Robjs/liver_kidney_atger_nestle/systems_clockdriven_tissuewide_genes.Robj", v=T)
  
  # jlambda <- 0.035  # liv only
  plda.out <- PenalizedLDA(M, M.labs, lambda = jlambda, K = 1, standardized = FALSE)
  
  # plot pretty
  vec.length <- sqrt(plda.out$discrim[, 1]^2)
  
  jsize.cutoff <- 0
  jsize.pairs.cut <- sapply(vec.length, function(jsize){
    if (jsize > jsize.cutoff){
      return(jsize)
    } else {
      return(0)
    }
  })
  
  labels <- names(plda.out$x)
  labels.cut <- mapply(function(jlab, jsize){
    if (jsize <= 0 & jlab != "SRF"){
      return("")
    } else {
      return(jlab)
    }
  }, labels, jsize.pairs.cut)
  
  dat.plot <- data.frame(discrim = plda.out$discrim[, 1],
                         motif = labels.cut,
                         motif.orig = labels,
                         vec.length = vec.length,
                         vec.length.cut = jsize.pairs.cut)
  dat.plot$discrim.floor <- sapply(dat.plot$discrim, function(d){
    if (d > 0){
      return("Systemic")
    } else if (d < 0){
      return("Clock")
    } else {
      return(NA)
    }
  })
  dat.plot <- dat.plot %>%
    arrange(discrim) %>%
    mutate(Index = seq(length(discrim)))
  dat.labs <- subset(dat.plot, vec.length.cut > 0)
  
  gene.plot <- data.frame(proj = plda.out$xproj, 
                          gene = rownames(plda.out$x),
                          jlabel = ifelse(plda.out$y == 1, "Clock", "Systems"))
  return(list(dat.plot = dat.plot, gene.plot = gene.plot))
}


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
save.image(file = "/data/shared/jake_data/tissue_specificity/PostPhDFiles/GR_2018_Primetime_Objects.Rdata")