# 2017-04-28
# Jake Yeung
# Can we identify rhythmic regulators by filtering liver DHSs then running MARA?

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

start.global <- Sys.time()

library(ggplot2)
library(reshape2)
library(dplyr)

source("scripts/functions/ListFunctions.R")
source("scripts/functions/HardcodedConstants.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")

GetTissSpecPeaks <- function(S.long, jgenes, distfilt, jcutoff, jcutoff.low, rhyth.tiss, flat.tiss){
  # get tiss spec peaks
  S.sub <- subset(S.long, gene %in% jgenes & dist < distfilt)
  jpeaks <- as.character(unique(S.sub)$peak)  # 192726 peaks for Liver genes within 50 kb away
  print(paste("number of peaks surrounding genes", length(jpeaks)))
  
  # take peaks with Liver signal greater than cutoff
  jtiss <- levels(S.sub$tissue)
  tiss.i <- which(jtiss %in% rhyth.tiss)
  others.i <- which(jtiss %in% flat.tiss)
  
  S.sub.tisspeaks <- S.sub %>%
    group_by(peak, gene) %>%  # tissue order as "Cere", "Heart", "Kidney", "Liver", "Lung", "Mus"
    filter(min(zscore[tiss.i]) > jcutoff & max(zscore[others.i]) < jcutoff.low) 

  return(S.sub.tisspeaks) 
}

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

# # Hard coded constants
# jweight <- 0  # take all liver DHSs to all genes in Liver_SV129
# distfilt <- 40000
# include.promoters <- FALSE
# all.genes <- FALSE
# use.sql <- TRUE
# jmod <- "Liver_SV129"
# # jmod <- "Liver_SV129,Liver_BmalKO"
# jmodstr <- gsub(",", "-", jmod)
# jcutoff <- 3
# jcutoff.low <- 0
# do.pairs <- TRUE

# Take from command line
args <- commandArgs(trailingOnly = TRUE)
print(args)
jweight <- check.type(args[[1]], checkfn = as.numeric)
distfilt <- check.type(args[[2]], checkfn = as.numeric)
all.genes <- FALSE
use.sql <- TRUE
jmod <- check.type(args[[3]], checkfn = as.character)
jcutoff <- check.type(args[[4]], checkfn = as.numeric)
jcutoff.low <- check.type(args[[5]], checkfn = as.numeric)
include.promoters <- check.type(args[[6]], checkfn = as.logical)
do.pairs <- check.type(args[[7]], checkfn = as.logical)
if (include.promoters){
  distfilt.min <- 1000  # if using include.promoters, need to exclude peaks that are too close to promoter of genes.
} else {
  distfilt.min <- 0
}
# include.promoters <- FALSE

if (jmod != "flat"){
  jmodstr <- gsub(",", "-", jmod)
} else {
  jmod <- ""
  jmodstr <- "flat"
}
jcutoffstr <- paste(jcutoff, jcutoff.low, sep = ".")
# suffix <- paste0(".weight.", jweight, ".promoters.", include.promoters, ".all_genes.", all.genes, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr)

suffix <- GetSuffix(jweight, use.sql, jmodstr, jcutoffstr, include.promoter = include.promoters, mindist = distfilt.min, do.pairs = do.pairs)
# suffix <- paste0(".weight.", jweight, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr)

# determine Rhyth tiss, Flat tiss programmatically
tiss <- c("Liver", "Kidney")
if (jmodstr != "flat"){
  rhyth.tiss <- strsplit(strsplit(jmod, ",")[[1]][[1]], split = "_")[[1]][[1]]
  if (!rhyth.tiss %in% tiss){
    stop(paste(jmod, rhyth.tiss, "must be Liver or Kidney"))
  }
  flat.tiss <- tiss[!tiss %in% rhyth.tiss]
} else {
  rhyth.tiss <- "Liver"
  flat.tiss <- "Kidney"
}
print(paste("Tissues:", rhyth.tiss, flat.tiss))

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)

liver.genes.all <- as.character(subset(fits.long.filt, model %in% jmod & weight >= jweight)$gene)

if (all.genes){
  liver.genes <- liver.genes.all
} else {
  if (!exists("S.long")) load("Robjs/S.long.multigene.filt.50000.Robj", v=T)
  S.sub.livpeaks <- GetTissSpecPeaks(S.long = S.long, jgenes = liver.genes.all, distfilt = distfilt, 
                                     jcutoff = jcutoff, jcutoff.low = jcutoff.low, rhyth.tiss = rhyth.tiss, flat.tiss = flat.tiss)
  liver.peaks <- unique(as.character(S.sub.livpeaks$peak))
  liver.genes <- unique(as.character(S.sub.livpeaks$gene))
  print(paste("N peaks:, ", length(liver.peaks)))
  print(paste("N genes:, ", length(liver.genes)))
}


# Print out liver peaks for analysis later --------------------------------

peaksgenes <- data.frame(peak = as.character(S.sub.livpeaks$peak),
                         gene = as.character(S.sub.livpeaks$gene), stringsAsFactors = FALSE)
peaksdir <- "/home/yeung/data/tissue_specificity/tissuepeaksgenes"
dir.create(peaksdir)
save(peaksgenes, file = file.path(peaksdir, paste0("liver.spec.peaks", suffix, ".Robj")))


# Get gene expression over time and genotypes -----------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)

# use "all" genes because filtering comes at sitecounts
dat.wtko <- subset(dat.wtko, gene %in% liver.genes.all)

dat.wtko$sampname <- paste(as.character(dat.wtko$tissue), sprintf("%02d", as.numeric(dat.wtko$time)), as.character(dat.wtko$geno), sep = "_") 
dat.mat <- dcast(dat.wtko, formula = "gene ~ sampname", value.var = "exprs")
dat.mat <- dplyr::rename(dat.mat, "Gene.ID" = gene)

do.center <- TRUE
if (do.center){
  # first column is gene name ignore it
  dat.mat[, -1] <- t(scale(t(dat.mat[, -1]), center = TRUE, scale = FALSE))
}

# Get sitecounts ----------------------------------------------------------

if (!use.sql){
  load("Robjs/liver_kidney_atger_nestle/N.long.all_genes.3wtmodules.bugfixed.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
  
  # N.long.filt <- subset(N.long.filt, gene %in% liver.genes & dist <= distfilt)
  N.long.filt <- subset(N.long.filt, gene %in% liver.genes & peak %in% liver.peaks & dist <= distfilt)
  
  print(length(unique(as.character(N.long.filt$gene))))
  print(length(unique(as.character(N.long.filt$peak))))
} else {
  inf <- "/home/shared/sql_dbs/closestbed_multiple_genes.genomewide.merged.motifindexed.sqlite3"
  motevo.tbl <- LoadDatabase(inf)
  
  print("Getting genes from database")
  start <- Sys.time()
  N.sub.lst <- expandingList()
  for (jgene in liver.genes){
    N.long.filt.query <- filter(motevo.tbl, gene == jgene)  # peaks are not indexed, so dont take them
    N.sub.tmp <- collect(N.long.filt.query, n = Inf)
    N.sub.lst$add(N.sub.tmp)
  }
  N.long.filt <- N.sub.lst$as.list()
  N.long.filt <- bind_rows(N.long.filt)
  
  rm(N.sub.tmp, N.sub.lst)  # worth it? 
  
  # filter peaks after querying database
  N.long.filt <- subset(N.long.filt, peak %in% liver.peaks & dist <= distfilt)
  print(paste("Collected", length(unique(as.character(N.long.filt$gene))), "genes and ", length(unique(as.character(N.long.filt$peak))), "peaks"))
  print(str(N.long.filt))
  print(Sys.time() - start)
  N.long.filt$motif <- make.names(N.long.filt$motif)
  
  if (do.pairs){
    source("scripts/functions/LdaFunctions.R")
    rhyth.motifs <- c("RORA.p2", "NFIL3.p2", "AR.p2", "bHLH_family.p2", "HSF1.2.p2", "SRF.p3", "NR3C1.p2")
    tissue.motifs <- c("ONECUT1.2.p2", "CUX2.p2", "FOXA2.p3")
    jsub <- subset(N.long.filt, motif %in% c(rhyth.motifs, tissue.motifs))
    jmat <- dcast(jsub, formula = "gene + dist + peak ~ motif", fill = 0, value.var = "sitecount", fun.aggregate = sum)
    jmat.meta <- subset(jmat, select = c(gene, dist, peak))
    jmat.exprs <- subset(jmat, select = -c(gene, dist, peak))
    jmat.rhyth <- jmat[, rhyth.motifs]
    jmat.tiss <- jmat[, tissue.motifs]
    # do an "OR" by summing
    # jmat.tiss <- data.frame(apply(jmat.tiss, 1, sum)); colnames(jmat.tiss) <- "LiverMotif.p2"
    jmat.cross <- CrossProductTwoSets(jmat.rhyth, jmat.tiss)
    # filter out 0s
    rows.keep <- apply(jmat.cross, 1, max) > 0
    jmat.crossmeta <- cbind(jmat.meta[rows.keep, ], jmat.cross[rows.keep, ])
    jmat.crossmeta$chromo <- "cross"
    jmat.crossmeta$start <- "cross"
    jmat.crossmeta$end <- "cross"
    jmat.crossmeta[jmat.crossmeta == 0] <- NA
    jmat.long <- melt(jmat.crossmeta, id.vars = c("chromo", "start", "end", "gene", "dist", "peak"), variable.name = "motif", value.name = "sitecount", na.rm = TRUE)
    N.long.filt <- rbind(N.long.filt, jmat.long)
    # N.long.filt <- rbind(subset(N.long.filt, !motif %in% c(rhyth.motifs, tissue.motifs)), jmat.long)
  }
}

if (include.promoters){
  # add promoter counts
  Npath <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
  N <- read.table(Npath, header=TRUE)
  N.prom.long <- melt(N, id.vars = "Gene.ID", variable.name = "motif", value.name = "sitecount")
  N.prom.long <- dplyr::rename(N.prom.long, "gene" = Gene.ID)
  N.prom.long$chromo <- NA
  N.prom.long$start <- NA
  N.prom.long$end <- NA
  N.prom.long$dist <- 0
  N.prom.long$peak <- "promoter"
 
  N.sum <- bind_rows(N.prom.long, subset(N.long.filt, dist >= distfilt.min)) %>%
    group_by(gene, motif) %>%
    summarise(sitecount = sum(sitecount))
  
} else {
  
  N.sum <- N.long.filt %>%
    group_by(gene, motif) %>%
    summarise(sitecount = sum(sitecount))
  
}
N.mat <- dcast(N.sum, formula = "gene ~ motif", value.var = "sitecount", fill = 0)
N.mat <- dplyr::rename(N.mat, "Gene.ID" = gene)

if (jmodstr == "flat"){
  stop("Stopping because for flat genes I only want to wake up the database")
}

# Add pairs to N  ---------------------------------------------------------

# if (do.pairs){
#   source("/home/yeung/projects/tissue-specificity/scripts/functions/LdaFunctions.R")
#   rhyth.motifs <- c("RORA.p2", "NFIL3.p2", "AR.p2", "bHLH_family.p2", "HSF1.2.p2", "SRF.p3", "NR3C1.p2")
#   tissue.motifs <- c("ONECUT1.2.p2", "CUX2.p2", "FOXA2.p3")
#   N.rhyth <- N.mat[, rhyth.motifs]
#   N.tiss <- N.mat[, tissue.motifs]
#   # do OR 
#   N.tiss <- data.frame(apply(N.tiss, 1, sum)); colnames(N.tiss) <- "LiverMotif.p2"
#   
#   N.rhythtiss <- CrossProductTwoSets(N.rhyth, N.tiss)
#   N.mat <- cbind(N.mat, N.rhythtiss)
# }


# Write E and N to output -------------------------------------------------


# write E
# E.subdir <- paste0("centered.", do.center, ".mod.", jmodstr, ".weightcutoff.", jweight)
E.subdir <- GetESubDir(do.center, jmodstr, jweight)
E.dir <- file.path("/home/yeung/data/tissue_specificity/gene_exprs_for_mara", E.subdir)
dir.create(E.dir, recursive = TRUE, showWarnings = FALSE)

# write Sample by Sample
tissues <- c("Liver", "Kidney")
genos <- c("SV129", "BmalKO")
sampnames.all <- colnames(dat.mat)  # include Gene.ID in your grep

for (tiss in tissues){
  for (geno in genos){
    cond <- paste(tiss, geno, sep = "_")
    sampnames.i <- grepl(paste0("Gene.ID|", tiss, "_[0-9]*_", geno), sampnames.all)
    E.out <- file.path(E.dir, paste0(cond, ".centered.", do.center, ".mat"))
    if (!file.exists(E.out)){
      write.table(dat.mat[, sampnames.i], E.out, append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } else {
      print("File exists, skipping writing of gene expression matrix")
    }
  }
}

# write N
N.dir <- paste0("/home/yeung/data/tissue_specificity/sitecount_matrices_for_mara")
N.out <- file.path(N.dir, paste0("sitecounts_swiss_regulon_promoters_only.weight", suffix, ".mat"))

dir.create(N.dir)

if (!file.exists(N.out)){
  write.table(N.mat, file = N.out, append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
} else {
  print("Sitecounts exists, skipping")
}


# Run MARA ----------------------------------------------------------------

# marascript <- "/home/yeung/projects/tissue-specificity/scripts/three_way_contingency_tables/run_run_filter_mara_pipeline.promoters.bygenelist.sh"
# requires E.dir (contains .Mat of exprs), N.mat (sitecounts, filtered), outmain (needs to be unique new folder)
marascript <- "/home/yeung/projects/sleep_deprivation/scripts/shellscripts/run_mara_sleep_deprivation.sh"
jmain <- "/home/yeung/data/tissue_specificity/mara_results"
dir.create(jmain)
outmain <- paste0(jmain, "/mara_outputs", suffix)
dir.create(outmain, showWarnings = FALSE)
outdir <- file.path(outmain, paste0("center.", do.center, suffix))
# do not create outdir, because MARA will not run on an already-created directory 

if(!dir.exists(outdir)){
  cmd <- paste("bash", marascript, N.out, outdir, E.dir)
  system(cmd)
  
  # consolidate
  combinescript <- "/home/yeung/projects/ridge-regression/run_scripts/combine_activities_and_plot.one_dir.sh"
  outdir.mara <- file.path(outdir, E.subdir)
  
  cmd2 <- paste("bash", combinescript, outdir.mara)
  system(cmd2)
}

print(Sys.time() - start.global)
