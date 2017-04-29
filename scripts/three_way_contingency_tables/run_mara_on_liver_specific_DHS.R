# 2017-04-28
# Jake Yeung
# Can we identify rhythmic regulators by filtering liver DHSs then running MARA?

rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)


# Load liver DHS ----------------------------------------------------------

jweight <- 0.8  # take all liver DHSs to all genes in Liver_SV129
jweight <- 0  # take all liver DHSs to all genes in Liver_SV129

flatampmax <- 0.1

promoters.only <- FALSE
all.genes <- FALSE
suffix <- paste0(".weight.", jweight, ".flatampmax.", flatampmax, ".promoters.", promoters.only, ".all_genes.", all.genes)

inf <- paste0("/home/yeung/projects/tissue-specificity/Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.200.weight.", jweight, ".MergePeaks.FALSE.nullmodel.JI.flatampmax.", flatampmax, ".withNmatallNmatfreqs.RemoveZeroCounts.Robj")
load(inf, v=T)

jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)

liver.genes.all <- as.character(subset(fits.long.filt, model == "Liver_SV129")$gene)

if (all.genes){
  liver.genes <- liver.genes.all
} else {
  liver.peaks <- unique(as.character(subset(N.mat.all, model %in% "rhyth")$peak))
  liver.genes <- unique(as.character(subset(N.mat.all, model %in% "rhyth")$gene))
  print(paste("N peaks:, ", length(liver.peaks)))
  print(paste("N genes:, ", length(liver.genes)))
}


# Get gene expression over time and genotypes -----------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)

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

# write to output
E.dir <- file.path("/home/yeung/data/tissue_specificity/gene_exprs_for_mara", paste0("centered.", do.center))
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


# Get sitecounts ----------------------------------------------------------

if (!promoters.only){
  distfilt <- 40000
  load("Robjs/liver_kidney_atger_nestle/N.long.all_genes.3wtmodules.bugfixed.Robj", v=T); N.long.filt <- N.long.livertwflat; rm(N.long.livertwflat)
  
  # N.long.filt <- subset(N.long.filt, gene %in% liver.genes & dist <= distfilt)
  N.long.filt <- subset(N.long.filt, gene %in% liver.genes & peak %in% liver.peaks & dist <= distfilt)
  
  print(length(unique(as.character(N.long.filt$gene))))
  print(length(unique(as.character(N.long.filt$peak))))
  
  N.sum <- N.long.filt %>%
    group_by(gene, motif) %>%
    summarise(sitecount = sum(sitecount))
  
  N.mat <- dcast(N.sum, formula = "gene ~ motif", value.var = "sitecount", fill = 0)
  N.mat <- dplyr::rename(N.mat, "Gene.ID" = gene)
} else {
  Npath <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
  N <- read.table(Npath, header=TRUE)
  N.mat <- subset(N, Gene.ID %in% liver.genes)
}

# write to output
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
cmd <- paste("bash", marascript, N.out, outdir, E.dir)
system(cmd)

# consolidate
combinescript <- "/home/yeung/projects/ridge-regression/run_scripts/combine_activities_and_plot.one_dir.sh"
outdir.mara <- file.path(outdir, "centered.TRUE")

cmd2 <- paste("bash", combinescript, outdir.mara)
system(cmd2)



