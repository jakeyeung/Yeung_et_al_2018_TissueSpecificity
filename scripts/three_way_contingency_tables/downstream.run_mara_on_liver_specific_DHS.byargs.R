# 2017-04-28
# Jake Yeung
# Can we identify rhythmic regulators by filtering liver DHSs then running MARA?

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")


library(ggplot2)
library(reshape2)
library(dplyr)


source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/CosSineFunctions.R")
source("scripts/functions/HalfLifeFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/ProteomicsFunctions.R")


check.type <- function(arg, checkfn){
  # force argument by checkfn, complain if NA
  arg.new <- checkfn(arg)
  if (!is.na(arg.new)){
    return(arg.new)
  } else {
    stop(paste("Arg:", arg, "did not pass check function"))
    return(NA)
  }
}


# Load liver DHS ----------------------------------------------------------

# # Hard coding
# jweight <- 0.8  # take all liver DHSs to all genes in Liver_SV129
# jweight <- 0  # take all liver DHSs to all genes in Liver_SV129
# do.center <- TRUE
# promoters.only <- FALSE
# all.genes <- FALSE
# use.sql <- TRUE
# jmod <- "Liver_SV129,Liver_BmalKO"
# jmod <- "Liver_SV129"
# jcutoff <- 3
# jcutoff.low <- 0

do.center <- TRUE

# Take from command line
args <- commandArgs(trailingOnly = TRUE)
print(args)
jweight <- check.type(args[[1]], checkfn = as.numeric)
distfilt <- check.type(args[[2]], checkfn = as.numeric)
# promoters.only <- check.type(args[[3]], checkfn = as.logical)
promoters.only <- FALSE
all.genes <- FALSE
use.sql <- TRUE
jmod <- check.type(args[[3]], checkfn = as.character)
jcutoff <- check.type(args[[4]], checkfn = as.numeric)
jcutoff.low <- check.type(args[[5]], checkfn = as.numeric)

jmodstr <- gsub(",", "-", jmod)
jtiss.nogeno <- strsplit(strsplit(jmod, ",")[[1]], "_")[[1]][[1]]
jcutoffstr <- paste(jcutoff, jcutoff.low, sep = ".")
# suffix <- paste0(".weight.", jweight, ".flatampmax.", flatampmax, ".promoters.", promoters.only, ".all_genes.", all.genes, ".sql.", use.sql)
# suffix <- paste0(".weight.", jweight, ".flatampmax.", flatampmax, ".promoters.", promoters.only, ".all_genes.", all.genes, ".sql.", use.sql, ".mod.", jmod)
# suffix <- paste0(".weight.", jweight, ".flatampmax.", flatampmax, ".promoters.", promoters.only, ".all_genes.", all.genes, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr)

# suffix <- paste0(".weight.", jweight, ".promoters.", promoters.only, ".all_genes.", all.genes, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr)
# suffix <- paste0(".weight.", jweight, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr)
suffix <- GetSuffix(jweight, use.sql, jmodstr, jcutoffstr)

# E.subdir <- paste0("centered.", do.center, ".mod.", jmodstr)
# E.subdir <- paste0("centered.", do.center, ".mod.", jmodstr, ".weightcutoff.", jweight)
E.subdir <- GetESubDir(do.center, jmodstr, jweight)


jmain <- "/home/yeung/data/tissue_specificity/mara_results"
outmain <- paste0(jmain, "/mara_outputs", suffix)
outdir <- file.path(outmain, paste0("center.", do.center, suffix))
maraoutdir <- file.path(outdir, E.subdir)

if (!dir.exists(maraoutdir)){
  print(paste(maraoutdir, "does not exist, exiting quietly"))
  quit(save = "no", status = 0, runLast = FALSE)
}

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)


# Get gene expression over time and genotypes -----------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)


# Load other stuff --------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
prot.long <- LoadProteomicsData()
prot.long <- subset(prot.long, geno == "WT" & tissue == "Liver")

# Load MARA output --------------------------------------------------------

# maraoutdir <- paste0("/home/yeung/data/tissue_specificity/mara_results/mara_outputs", suffix, "/center.TRUE", suffix, "/centered.TRUE")


act.s <- LoadActivitiesLong(indir = maraoutdir, shorten.motif.name = TRUE, make.cnames = FALSE)
act.s$sampname <- act.s$tissue
act.s$tissue <- as.character(sapply(as.character(act.s$sampname), function(s) strsplit(s, "_")[[1]][[1]]))
act.s$time <- as.numeric(sapply(as.character(act.s$sampname), function(s) strsplit(s, "_")[[1]][[2]]))
act.s$geno <- as.character(sapply(as.character(act.s$sampname), function(s) strsplit(s, "_")[[1]][c(-1, -2)]))
act.s$tissue <- paste(act.s$tissue, act.s$geno, sep = "_")
act.s$tissue <- factor(act.s$tissue, levels = c("Liver_SV129", "Liver_BmalKO", "Kidney_SV129", "Kidney_BmalKO"))
act.s$experiment <- "rnaseq"
act.s$sampname <- NULL

fourier.scale <- 4
zscore.min <- 1.25
omega <- 2 * pi / 24
hr.shift <- 3
mrna.hl <- log(2) / (omega / tan(omega * hr.shift))  # convert hour shift to half-life 
comp <- 1

act.s.complex <- ProjectWithZscore(act.s, omega, n = fourier.scale)
sig.motifs <- unique(as.character(subset(act.s.complex, zscore > zscore.min)$gene))
s.act <- SvdOnComplex(subset(act.s.complex, gene %in% sig.motifs), value.var = "exprs.transformed")

eigens.act <- GetEigens(s.act,
                        eigenval = FALSE,
                        period = 24, 
                        comp = comp, 
                        adj.mag = TRUE, 
                        constant.amp = 6, 
                        label.n = 20, 
                        jtitle = "", 
                        peak.to.trough = TRUE, 
                        dot.col = "black", 
                        dotsize = 2, 
                        dotshape = 18, 
                        label.gene = c("ELF1.2.4"),
                        half.life = mrna.hl)

eigens.act.fancy.LivWTKO <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, 
                                      eigenval = FALSE,
                                      constant.amp = 5, 
                                      label.n = Inf, jtitle = "", 
                                      peak.to.trough = TRUE, 
                                      dot.col = "black", 
                                      # dot.col = col.hash, 
                                      dotsize = 6, 
                                      dotshape = 18,
                                      disable.text = FALSE, 
                                      add.arrow = TRUE,
                                      disable.repel = TRUE,
                                      half.life = mrna.hl)



# also plot TF activity adjusted by half-life
# s.ampphase <- GetAmpPhaseFromActivities(act.s, mrna.hl, jtiss = jtiss, jgeno = "SV129")
# just shift by 3 hours
act.s.shift <- act.s
act.s.shift$time <- act.s$time - hr.shift
act.s.shift$time <- sapply(act.s.shift$time, function(x) ifelse(x < 0, x + 48, x))

jmotif <- "CEBPA.B_DDIT3"
jmotif <- "ESR1"
jmotif <- "RORA"
jmotif <- "DBP"
print(PlotActivitiesWithSE(subset(act.s.shift, gene == jmotif), jtitle = jmotif) + theme_bw())





jtiss <- strsplit(jmod, ",")[[1]][[1]]
print(jtiss)
liv.rhyth <- as.character(subset(fits.bytiss, tissue == jtiss & amp > 0.1 & pval < 1e-2)$gene)
if (length(liv.rhyth) == 0){
  warning("No rhythmic genes!")
}
tfs <- GetTFs(get.mat.only = TRUE)

# Print genes that match model
jmotifs <- names(head(eigens.act$eigensamp[order(abs(eigens.act$eigensamp), decreasing = TRUE)], n = 20))


# Plot things


pdf(paste0("/home/yeung/projects/tissue-specificity/plots/mara_liver_kidney_modules_on_liverDHS/plots", suffix, ".pdf"))

# print(eigens.act$v.plot)
print(eigens.act$u.plot)

# print(eigens.act.fancy.LivWTKO$v.plot)
print(eigens.act.fancy.LivWTKO$u.plot)



for (jmotif in jmotifs){
  
  print(PlotActivitiesWithSE(subset(act.s.shift, gene == jmotif), jtitle = jmotif) + theme_bw() + theme(aspect.ratio = 1))
  
  genes.all <- unlist(sapply(jmotif, GetGenesFromMotifs, tfs))
  genes.that.fit <- genes.all[which(genes.all %in% liv.rhyth)]
  if (jmotif == "RORA"){
    genes.that.fit <- c(genes.that.fit, "Nr1d1")
  }
  print(paste(jmotif, ": Rhyth genes:", paste(genes.that.fit, collapse = ",")))
  if (length(genes.that.fit) > 0){
    for (gene.hit in genes.that.fit){
      if (jmod %in% c("Liver_SV129,Liver_BmalKO", "Liver_SV129", "Liver_BmalKO")){
        gene.prot <- gene.hit
        jprot.long <- prot.long
      } else {
        gene.prot <- ""
        jprot.long <- NA
      }
      print(PlotGeneTissuesWTKO(subset(dat.wtko, gene == gene.hit), jtitle = gene.hit))
      print(PlotmRNAActivityProtein(dat.wtko, act.s.shift, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss.nogeno, dotsize = 3, themesize = 22, single.day = TRUE) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, act.s, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, act.s.shift, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, s.ampphase, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22, act.in.sine.cos = TRUE) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, s.ampphase, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = jtiss, dotsize = 3, themesize = 22, act.in.sine.cos = TRUE, single.day = TRUE) + theme(strip.text = element_blank()))
      # print(PlotmRNAActivityProtein(dat.wtko, act.s, gene.dat = gene.hit, prot.long = jprot.long, gene.act = jmotif, gene.prot = gene.prot, jtiss = "both", dotsize = 2, themesize = 14) + theme(strip.text = element_blank()))
    }
  }
}
dev.off()

# 
# # Old ---------------------------------------------------------------------
# 
# 
# mara <- LoadMaraOutput(maraoutdir)
# 
# zscores <- mara$zscores
# act.long <- mara$act.long
# act.long$samp <- NULL
# 
# act.long$tissue <- sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[1]])
# act.long$time <- as.numeric(sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[2]]))
# act.long$geno <- as.character(sapply(as.character(act.long$sample), function(s) strsplit(s, "_")[[1]][[3]]))
# 
# PlotGeneTissuesWTKO(subset(act.long, gene == "RORA.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "bHLH_family.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "TFAP4.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "HNF1A.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "FOSL2.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "AR.p2"))
# PlotGeneTissuesWTKO(subset(act.long, gene == "FOXA2.p3"))
# 
# 
