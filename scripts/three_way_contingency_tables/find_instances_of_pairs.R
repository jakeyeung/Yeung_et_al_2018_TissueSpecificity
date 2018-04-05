# 2017-05-05
# Jake Yeung
# find pairs: first get cutoffs for each TF, and then use cutoff to find "hits"
# plot "hits" on genome-browser.

rm(list=ls())

library(hash)
library(dplyr)
library(ggplot2)
library(plotrix)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/ListFunctions.R")
source("scripts/functions/HardcodedConstants.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/HandleMotifNames.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/CosSineFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LiverKidneyFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/ProteomicsFunctions.R")
source("scripts/functions/FisherTestSitecounts.R")
source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")



# Functions ---------------------------------------------------------------

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


# Constants ---------------------------------------------------------------


do.center <- TRUE

distfilt <- 40000
jweight <- 0.8
use.sql <- TRUE
jmod <- "Liver_SV129,Liver_BmalKO"
jmod <- "Liver_SV129"
# args <- commandArgs(trailingOnly = TRUE)
# jcutoff <- as.numeric(args[[1]])
jcutoff <- 3.0
jcutoff.low <- 0
incl.promoters <- FALSE

rhyth.tiss <- "Liver"
flat.tiss <- "Kidney"

jcutoffstr <- paste(jcutoff, jcutoff.low, sep = ".")
jmodstr <- gsub(",", "-", jmod)
jtissgeno <- strsplit(jmod, ",")[[1]][[1]]
jtiss.nogeno <- strsplit(strsplit(jmod, ",")[[1]], "_")[[1]][[1]]

suffix <- GetSuffix(jweight, use.sql, jmodstr, jcutoffstr, incl.promoters)
E.subdir <- GetESubDir(do.center, jmodstr, jweight)


# Load data ---------------------------------------------------------------

# get flat genes and liver genes
jmeth <- "g=1001"
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)

liver.genes.all <- as.character(subset(fits.long.filt, model %in% jmod & weight >= jweight)$gene)
flat.genes.all <- as.character(subset(fits.long.filt, model %in% "" & weight >= jweight)$gene)

load("Robjs/S.long.multigene.filt.50000.Robj", v=T)

# get liver peaks at liver genes
S.sub.livpeaks <- GetTissSpecPeaks(S.long = S.long, jgenes = liver.genes.all, distfilt = distfilt,
                                   jcutoff = jcutoff, jcutoff.low = jcutoff.low, rhyth.tiss = rhyth.tiss, flat.tiss = flat.tiss)
liver.peaks <- unique(as.character(S.sub.livpeaks$peak))
liver.genes <- unique(as.character(S.sub.livpeaks$gene))
print(paste("N peaks:, ", length(liver.peaks)))
print(paste("N genes:, ", length(liver.genes)))

# get flat
S.sub.flat <- GetTissSpecPeaks(S.long = S.long, jgenes = flat.genes.all, distfilt = distfilt,
                                   jcutoff = jcutoff, jcutoff.low = jcutoff.low, rhyth.tiss = rhyth.tiss, flat.tiss = flat.tiss)
flat.peaks <- unique(as.character(S.sub.flat$peak))
flat.genes <- unique(as.character(S.sub.flat$gene))

# get kidney peaks at liver genes
S.sub.kidpeaks <- GetTissSpecPeaks(S.long = S.long, jgenes = liver.genes.all, distfilt = distfilt,
                                   jcutoff = jcutoff * 0.5, jcutoff.low = jcutoff.low, rhyth.tiss = flat.tiss, flat.tiss = rhyth.tiss)
kid.peaks <- unique(as.character(S.sub.kidpeaks$peak))
kid.genes <- unique(as.character(S.sub.kidpeaks$gene))

# # if you load it from Robj
# maindir <- "/home/yeung/data/tissue_specificity/tissuepeaksgenes"
# prefix <- "liver.spec.peaks"
# inf1 <- file.path(maindir, paste0(prefix, suffix, ".Robj"))
# load(inf1, verbose = T)
# liver.peaks <- unique(as.character(peaksgenes$peak))
# liver.genes <- unique(as.character(peaksgenes$gene))
# # Get flat peaks
# suffix.flat <- GetSuffix(jweight, use.sql, "flat", jcutoffstr, incl.promoters)
# inf.flat <- file.path(maindir, paste0(prefix, suffix.flat, ".Robj"))
# load(inf.flat, v=T)
# flat.peaks <- unique(as.character(peaksgenes$peak))
# flat.genes <- unique(as.character(peaksgenes$gene))


# Load database -----------------------------------------------------------

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
N.long.liverpeaks <- subset(N.long.filt, peak %in% liver.peaks & dist <= distfilt)
N.long.kidpeaks <- subset(N.long.filt, peak %in% kid.peaks & dist <= distfilt)
N.long.liverpeaks$motif <- sapply(N.long.liverpeaks$motif, function(m) make.names(RemoveP2Name(m)))
N.long.kidpeaks$motif <- sapply(N.long.kidpeaks$motif, function(m) make.names(RemoveP2Name(m)))

print(paste("Collected liver peaks:", length(unique(as.character(N.long.liverpeaks$gene))), "genes and ", length(unique(as.character(N.long.liverpeaks$peak))), "peaks"))
print(paste("Collected kidney peaks:", length(unique(as.character(N.long.kidpeaks$gene))), "genes and ", length(unique(as.character(N.long.kidpeaks$peak))), "peaks"))
print(Sys.time() - start)



print("Getting flat genes from database")
start <- Sys.time()
N.sub.lst <- expandingList()
for (jgene in flat.genes){
  N.long.filt.query <- filter(motevo.tbl, gene == jgene)  # peaks are not indexed, so dont take them
  N.sub.tmp <- collect(N.long.filt.query, n = Inf)
  N.sub.lst$add(N.sub.tmp)
}
N.long.flat <- N.sub.lst$as.list()
N.long.flat <- bind_rows(N.long.flat)
N.long.flat <- subset(N.long.flat, peak %in% flat.peaks & dist <= distfilt)
# can do motif after filtering because you only filter once
N.long.flat$motif <- sapply(N.long.flat$motif, function(m) make.names(RemoveP2Name(m)))

# save(N.long.filt, N.long.flat, file="Robjs/N.long.flat.from.sql.Robj")
print(Sys.time() - start)

# label models
N.long.liverpeaks$model <- "Liver"
N.long.kidpeaks$model <- "Kidney"
N.long.flat$model <- "Flat"

N.merged <- bind_rows(N.long.liverpeaks, N.long.kidpeaks, N.long.flat)

N.merged <- N.merged %>%
  group_by(motif, gene, peak, model) %>%
  summarise(sitecount = sum(sitecount))

# fill in missing peaks with 0s
N.mat <- dcast(N.merged, formula = "gene + peak + model ~ motif", value.var = "sitecount", fill = 0)
N.merged <- melt(N.mat, id.vars = c("gene", "peak", "model"), variable.name = "motif", value.name = "sitecount")

# Find cutoff for rhythmic factors

pdf(paste0("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/prob_cutoffs.", suffix, ".pdf"))
cutoffs <- seq(from = 0, to = 1, by = 0.1)
N.rhythfactors <- RunFisherOnPromoters(N.merged, foreground.models = "Liver", background.models = "Flat", cutoffs = cutoffs, return.full.df = TRUE)
# Cutoff at ~0.25
jmotif <- "RORA"
ggplot(subset(N.rhythfactors, motif == jmotif), aes(x = cutoff, y = -log10(p.value))) + geom_point() + geom_line() + ggtitle(jmotif)

# Find cutoff for liver factors
N.liverfactors <- RunFisherOnPromoters(N.merged, foreground.models = "Liver", background.models = "Kidney", cutoffs = cutoffs, return.full.df = TRUE)
jmotif <- "FOXA2"
ggplot(subset(N.liverfactors, motif == jmotif), aes(x = cutoff, y = -log10(p.value))) + geom_point() + geom_line() + ggtitle(jmotif)
jmotif <- "CUX2"
ggplot(subset(N.liverfactors, motif == jmotif), aes(x = cutoff, y = -log10(p.value))) + geom_point() + geom_line() + ggtitle(jmotif)
jmotif <- "ONECUT1.2"
ggplot(subset(N.liverfactors, motif == jmotif), aes(x = cutoff, y = -log10(p.value))) + geom_point() + geom_line() + ggtitle(jmotif)
jmotif <- "RORA"
ggplot(subset(N.liverfactors, motif == jmotif), aes(x = cutoff, y = -log10(p.value))) + geom_point() + geom_line() + ggtitle(jmotif)
jmotif <- "bHLH_family"
ggplot(subset(N.liverfactors, motif == jmotif), aes(x = cutoff, y = -log10(p.value))) + geom_point() + geom_line() + ggtitle(jmotif)
ggplot(subset(N.liverfactors, motif %in% c("FOXA2", "CUX2", "ONECUT1.2", "RORA")), aes(x = cutoff, y = -log10(p.value), colour = motif)) + geom_point() + geom_line() + ggtitle(jmotif)

N.toplivers <- RunFisherOnPromoters(N.merged, foreground.models = "Liver", background.models = "Kidney", cutoffs = c(0.4, 0.5, 0.6), return.full.df = FALSE)
dev.off()


# Find genes associated with RORA, FOXA2/ONECUT/CUX2 ----------------------

clock.motif <- "RORA"
jmotifs <- c("FOXA2", "CUX2", "ONECUT1.2", clock.motif)
N.sub <- subset(N.merged, motif %in% jmotifs & model == "Liver")

# ror.cutoff <- 0.5
ror.cutoff <- 0.8
# foxa2.cutoff <- 0.32
foxa2.cutoff <- 0.8
onecut.cutoff <- 0.8
cux2.cutoff <- 0.8
jcutoffs <- c(foxa2.cutoff, cux2.cutoff, onecut.cutoff, ror.cutoff)
cutoff.df <- data.frame(motif = jmotifs, 
                        cutoff = jcutoffs)

N.sub$above.cutoff <- mapply(function(jmotif, jsitecount) ifelse(jsitecount > subset(cutoff.df, motif == jmotif)$cutoff, TRUE, FALSE), as.character(N.sub$motif), N.sub$sitecount)

# filter TRUEs
N.sub.filt <- N.sub %>%
  group_by(gene, peak, model) %>%
  filter(above.cutoff == TRUE)

N.sub.filt$motif <- as.character(N.sub.filt$motif)
N.sub.filt$motif <- ifelse(N.sub.filt$motif == clock.motif, "Clock", "Tissue")

# Mark genes as containing ROR or Tissue
N.sub.filt.sum <- N.sub.filt %>%
  group_by(peak, gene) %>%  # or by gene only??
  # group_by(gene) %>%  # or by gene only??
  summarise(motif.hits = paste(unique(motif), collapse=",")) %>%
  group_by(gene) %>% 
  arrange(desc(motif.hits)) %>%
  filter(motif.hits != "Tissue") %>%
  filter(motif.hits == motif.hits[[1]])
  

# Plot polar plot, label genes
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.bytiss <- subset(fits.bytiss, tissue == "Liver_SV129" & gene %in% liver.genes)

motif.hits <- hash(N.sub.filt.sum$gene, N.sub.filt.sum$motif.hits)

# add to fits.bytiss as color
fits.bytiss$motif.hits <- sapply(as.character(fits.bytiss$gene), function(g){
  m.hits <- motif.hits[[g]]
  if (is.null(m.hits)){
    return(NA)
  } else {
    return(m.hits)
  }
})

# plot grays before colors by rearranging the dataframe
fits.bytiss <- fits.bytiss %>% arrange(desc(motif.hits)) %>% arrange(-row_number())

# find the Clock,Tissue or Tissue,Clock
tissclock.str <- unique(N.sub.filt.sum$motif.hits)[grepl(",", unique(N.sub.filt.sum$motif.hits))]

col.hash <- sapply(fits.bytiss$motif.hits, function(m){
  if (is.na(m)){
    jcol <- "gray85"
  } else if (m == "Tissue"){
    jcol <- "gray65"
  } else if (m == "Clock"){
    jcol <- "red"
  } else if (m == tissclock.str){
    jcol <- "darkred"
  }
})

fits.bytiss$label <- mapply(function(m, g) ifelse(m %in% c(tissclock.str, "Clock"), g, NA), fits.bytiss$motif.hits, as.character(fits.bytiss$gene))

plot.complex.nocol <- PlotComplex2(complex(modulus = fits.bytiss$amp * 2, argument = fits.bytiss$phase * 2 * pi / 24),
                                   ampscale = 1,
                                   labels = fits.bytiss$label,
                                   omega = 2 * pi / 24,
                                   add.arrow = FALSE,
                                   disable.repel = FALSE,
                                   constant.amp = 4,
                                   disable.text = FALSE,
                                   dot.col = col.hash,
                                   dotsize = 3,
                                   title = paste(clock.motif, suffix))

pdf(paste0("/home/yeung/projects/tissue-specificity/plots/liver_dhs_enrichment_plots/plot_complex_colored.", suffix, ".pdf"))
print(plot.complex.nocol)
dev.off()

# 
# # test on single motif
# jsub <- subset(N.merged.filled, motif == "RORA")
# 
# best.cutoff <- 0.3
# test <- FisherTestSitecounts(jsub, best.cutoff, show.table = TRUE)
# N.best <- RunFisherOnPromoters(N.merged.filled, foreground.models = "liver", background.models = "flat", cutoffs = best.cutoff, return.full.df = FALSE)
# 

# # Do enrichment analysis by hypergeometric test ---------------------------
# 
# N.gene <- N.long.filt %>%
#   group_by(gene, motif) %>%
#   summarise(sitecount = sum(sitecount))
#   
#   
# cutoffs <- seq(from = 0.4, to = 0.8, by = 0.1)
# N.sub <- RunFisherOnPromoters(N.long, foreground.models = models.tw, background.models = models.flat, cutoffs = cutoffs)
# 
# 
# # Load matrces ------------------------------------------------------------
# 
# load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T); dat.wtko <- dat.long; rm(dat.long)
# dat.wtko <- StaggeredTimepointsLivKid(dat.wtko)
# 
# load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
# fits.bytiss <- subset(fits.bytiss, tissue == "Liver_SV129" & gene != "")
# 
# prot.long <- LoadProteomicsData()
# prot.long <- subset(prot.long, geno == "WT" & tissue == "Liver")



