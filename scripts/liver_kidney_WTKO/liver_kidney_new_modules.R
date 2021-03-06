# 2016-07-26
# Jake Yeung
# liver_kidney_new_modules.R
# Find some exotic modules: maybe use that to find Kidney factors?


args.to.numeric <- function(arg, default = FALSE){
  # take command line argument and convert to numeric, if blank
  # then return default
  jarg <- as.numeric(arg)
  if (is.na(jarg)){
    return(default)
  } else {
    return(jarg)
  }
}


setwd("/home/yeung/projects/tissue-specificity")
jmeth <- "g=1001"
args <- commandArgs(trailingOnly = TRUE)
lambda <- args.to.numeric(args[1])

library(dplyr)
library(ggplot2)
library(hash)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/GetTFs.R")



# Functions ---------------------------------------------------------------



# Load --------------------------------------------------------------------


# load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
dat.orig <- dat.long

dat.long <- CollapseTissueGeno(dat.long)
dat.long <- StaggeredTimepointsLivKid(dat.long)

# filter NA changes
dat.long <- subset(dat.long, !is.na(gene))


# # # change Ciart to Gm129
dat.long$gene <- as.character(dat.long$gene)
dat.long$gene[which(dat.long$gene == "Ciart")] <- "Gm129"
print(head(dat.long))

# Annotate fits
fits.long.filt <- subset(fits.long.filt, method == jmeth)

# filter out genes
filt.genes <- c("")
fits.long.filt <- subset(fits.long.filt, !gene %in% filt.genes)
fits.long.filt$n.params <- sapply(fits.long.filt$model, function(m) return(length(strsplit(as.character(m), ";")[[1]])))
fits.long.filt$n.rhyth <- sapply(fits.long.filt$model, GetNrhythFromModel)
fits.long.filt$amp.avg <- mapply(GetAvgAmpFromParams, fits.long.filt$param.list, fits.long.filt$model)

fits.long.filt$gene <- as.character(fits.long.filt$gene)
fits.long.filt$gene[which(fits.long.filt$gene == "Ciart")] <- "Gm129"


# Init params -------------------------------------------------------------


suffix <- "globalLambda"  # playing
suffix <- ""  # normal

# lambda <- 0.0242444713967604
# lambda <- 0.0479929486250763  # tissue-wide clock-driven module
# lambda <- FALSE
# jmodel <- c("Liver_SV129,Kidney_SV129")
if (suffix != ""){
  suffix <- paste0(".", suffix)
}
# jmodel <- c("Liver_SV129,Kidney_SV129")
# jmodel <- c("Liver_SV129,Kidney_SV129,Liver_BmalKO,Kidney_BmalKO", "Liver_SV129,Liver_BmalKO;Kidney_SV129,Kidney_BmalKO")
# jmodel <- c("Kidney_SV129")
# jmodel <- c("Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO", "Liver_SV129,Kidney_SV129")
jmodel <- c("Liver_SV129;Liver_BmalKO")

if (jmodel == "Kidney_SV129"){
  # remove some outliers
  outliers <- c("Cox8b", "Ucp1", "Cidea", "Flg", "Nr4a2")
  fits.long.filt <- subset(fits.long.filt, !gene %in% outliers)
  dat.long <- subset(dat.long, !gene %in% outliers)
}

# # optionally using min rhyth
# min.rhyth <- 4
# jmodel <- as.character(unique(subset(fits.long.filt, n.rhyth >= min.rhyth)$model))
# # exclude models that may be problematic
# mod.excl <- "Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO"
# jmodel <- jmodel[which(jmodel != mod.excl)]


# Filter to common genes --------------------------------------------------

genes.keep <- unique(as.character(fits.long.filt$gene))

dat.long <- subset(dat.long, gene %in% genes.keep)

# Project to Frequency ----------------------------------------------------


omega <- 2 * pi / 24
dat.freq <- dat.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))

s <- SvdOnComplex(dat.freq, value.var = "exprs.transformed")

for (i in seq(1)){
  eigens <- GetEigens(s, period = 24, comp = i, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
}

# All periods -------------------------------------------------------------

periods <- rep(48, 6) / seq(1, 6)  # 48/1, 48/2 ... 48/12
loadfile <- "Robjs/liver_kidney_atger_nestle/dat.complex.all_T.Robj"
if (file.exists(loadfile)){
  load(loadfile)
} else {
  
  library(parallel)
  dat.complexes <- mclapply(periods, function(period, dat.long){
    omega <- 2 * pi / period
    dat.tmp <- dat.long %>%
      group_by(gene, tissue) %>%
      do(ProjectToFrequency2(., omega, add.tissue=TRUE))
    dat.tmp$period <- period
    return(dat.tmp)
  }, dat.long = dat.long, mc.cores = length(periods))
  
  
  dat.complex.all_T <- do.call(rbind, dat.complexes)
  outfcomp <- "Robjs/liver_kidney_atger_nestle/dat.complex.all_T.bugfixed.Robj"
  if (!file.exists(outfcomp)) save(dat.complex.all_T, file = outfcomp)
  rm(dat.complexes)
}
outffreq <- "Robjs/liver_kidney_atger_nestle/dat.freq.bugfixed.Robj"
if (!file.exists(outffreq)) save(dat.freq, file = outffreq)





i <- 1

fits.count <- subset(fits.long.filt, method == jmeth & model != "") %>% group_by(model) %>% summarise(model.count = length(model))
fits.count <- fits.count[order(fits.count$model.count, decreasing = TRUE), ]
fits.count <- subset(fits.count, model.count > 190)  # no underdetermination
jmodel.lst <- as.character(fits.count$model)


# for systems-driven module, all modules with >= 3 rhythmic conditions

jmodel.lst <- list(jmodel)
# jmodel.lst <- c("all")
for (jmodel in jmodel.lst){
  print(paste("Running model:"))
  print(jmodel)
  
  if (length(jmodel) <= 4){
    jmod <- paste(jmodel, collapse = "-")
    jmod <- gsub(";", "\\.", jmod)
    jmod <- paste0(jmod, ".", jmeth)
  } else {
    print("Shortening name to:")
    jmod <- paste0("many_modules_minrhyth.", min.rhyth, ".", jmeth)
  }
  
  if (jmodel != "all"){
    genes.tw <- as.character(subset(fits.long.filt, method == jmeth & model %in% jmodel)$gene)
  } else {
    genes.tw <- as.character(fits.long.filt$gene)
  }
  s <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = i, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE)
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
  
  
  # Write sitecounts --------------------------------------------------------
  
  Npath <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
  Npath.base <- paste0("/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/promoters_filtered_by_gene_liver_kidney_novel_modules", suffix)
  dir.create(Npath.base)
  Npath.out <- file.path(Npath.base, paste0("sitecount_matrix_geneids.", jmod, ".mat"))
  # Npath.out <- paste0("/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/promoters_filtered_by_gene_liver_kidney_novel_modules.", suffix, "/sitecount_matrix_geneids.", jmod, ".mat")
  N <- read.table(Npath, header=TRUE)
  
  N.sub <- subset(N, Gene.ID %in% genes.tw)
  print(paste("Writing", nrow(N.sub), "promoters to file:", Npath.out))
  
  if (!file.exists(Npath.out)) write.table(N.sub, file = Npath.out, append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  
  # Run MARA ----------------------------------------------------------------
  
  run.mara <- TRUE
  marascript <- "/home/yeung/projects/tissue-specificity/scripts/liver_kidney_WTKO/run_run_filter_mara_pipeline.promoters.bygenelist.singlemat.sh"
  nmat <- Npath.out
  
  outbase <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney", suffix)
  dir.create(outbase)
  if (!is.numeric(lambda)){
    outmain <- file.path(outbase, paste0("promoters.", jmod))
  } else {
    outmain <- file.path(outbase, paste0("promoters.", jmod, ".lambda.", lambda))
  }
  # outmain <- paste0("/home/yeung/projects/tissue-specificity/results/MARA.liver_kidney/promoters.", jmod)
  
  exprs.dir <- "/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/atger_with_kidney.bugfixed"
  cmd <- paste("bash", marascript, nmat, outmain, exprs.dir, lambda)
  if (!dir.exists(outmain)){
    print(cmd)
    system(cmd)
  }
  
  # Analyze output ----------------------------------------------------------
  
  indir <- file.path(outmain, "atger_with_kidney.bugfixed")
  source("scripts/functions/LoadActivitiesLong.R")
  act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = FALSE)
  # print(PlotActivitiesWithSE.rnaseq(subset(act.long, gene == jmotif & experiment == "rnaseq")) + theme_bw(24))
  # rename motifs based on the motifs with non-zero entries
  
  omega <- 2 * pi / 24
  act.complex <- act.long %>%
    group_by(gene, tissue) %>%
    do(ProjectToFrequency2(., omega, add.tissue=TRUE))
  
  s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")
  
  
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  # jtitle <- gsub(pattern = "\\.", replacement = "\n", basename(indirmain))
  
  max.labs <- 20
  jtitle <- ""
  for (comp in seq(1)){
    eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = max.labs, jtitle = jtitle, peak.to.trough = TRUE)
    print(eigens.act$u.plot + ggtitle(jmod))
    # multiplot(eigens.act$u.plot, eigens.act$v.plot, cols = 2)
  }
  
  tfs <- GetTFs(get.mat.only = TRUE)
  tfs$V2 <- as.character(tfs$V2)
  
  jmotifs <- names(head(eigens.act$eigensamp[order(abs(eigens.act$eigensamp), decreasing = TRUE)], n = max.labs))
  # take top and plot to output
  plotdir <- "plots/mara_liver_kidney_modules"
  dir.create(plotdir, showWarnings = FALSE)
  
  # BEGIN PLOT
  if (!is.numeric(lambda)){
    pdf(file.path(plotdir, paste0("motif_activity.", jmod, ".pdf")))
  } else {
    pdf(file.path(plotdir, paste0("motif_activity.", jmod, ".lambda.", lambda, ".pdf")))
  }
  
  multiplot(eigens$u.plot, eigens$v.plot, layout = jlayout)
  print(eigens.act$u.plot + ggtitle(jmod))
  
  for (jmotif in jmotifs){
    print(PlotActivitiesWithSE.rnaseq(subset(act.long, gene == jmotif & experiment == "rnaseq")) + theme_bw(24))
    genes <- GetGenesFromMotifs(jmotif, tfs)
    # genes <- strsplit(tfs[jmotif, ], ",")[[1]]
    # if (all(is.na(genes))){
    #   # try grep
    #   jgrep <- strsplit(jmotif, "\\.")[[1]][[1]]
    #   genes <- strsplit(tfs[grepl(jgrep, rownames(tfs)), ], ",")[[1]]
    # }
    print(paste(jmotif, ": ", length(genes), " genes"))
    for (g in genes){
      jsub <- subset(dat.orig, gene == g)
      if (nrow(jsub) > 0) print(PlotGeneTissuesWTKO(jsub, jtitle = g))
    }
  }
  
  # Print genes that match model
  genes.all <- unlist(sapply(jmotifs, GetGenesFromMotifs, tfs))
  genes.that.fit <- as.character(subset(fits.long.filt, gene %in% genes.all & model %in% jmodel & method == jmeth)$gene)
  if (length(genes.that.fit) > 0){
    for (gene.hit in genes.that.fit){
      print(gene.hit)
      print(PlotGeneTissuesWTKO(subset(dat.orig, gene == gene.hit), jtitle = gene.hit))
    }
  }
  
  dev.off()
}

