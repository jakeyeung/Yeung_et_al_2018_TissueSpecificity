# 2017-05-12
# Jake Yeung
# Hic1 and RORE do they always show up together?

rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(hash)

setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/DhsFunctions.R")
source("scripts/functions/GetClockGenes.R")
source("scripts/functions/ListFunctions.R")
source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")

FillWithZeros <- function(N.sub, jform = "gene + peak ~ motif", jval = "sitecount", idname = c("gene", "peak"), varname = "motif", valname = "sitecount"){
  jsub <- dcast(N.sub, formula = jform, value.var = jval, fill = 0)
  jsub <- melt(jsub, id.vars = idname, variable.name = varname, value.name = valname)
  return(jsub)
}

# Load stuff --------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)  # hogenesch
# model selection with g1000
load("Robjs/nconds_g1000_11_tissues/fits_long.11_tiss_3_max.g1000.bestmodel.filteramp.0.15.Robj", v=T)  # use this, looks same as fits.best should be OK?
load("Robjs/N.long.promoters_500.Robj", v=T)
load("Robjs/S.long.multigene.filt.50000.Robj", v=T)

jgenes <- GetClockGenes()
jgenes <- as.character(subset(fits.long, n.rhyth >= 8)$gene)
distfilt <- 40000
jcutoff <- 3
jcutoff.low <- Inf
rhyth.tiss <- "Liver"
flat.tiss <- "Kidney"

S.sub.peaks <- GetTissSpecPeaks(S.long = S.long, jgenes = jgenes, distfilt = distfilt,
                                   jcutoff = jcutoff, jcutoff.low = jcutoff.low, rhyth.tiss = rhyth.tiss, flat.tiss = flat.tiss)
jpeaks <- unique(as.character(S.sub.peaks$peak))

inf <- "/home/shared/sql_dbs/closestbed_multiple_genes.genomewide.merged.motifindexed.sqlite3"
motevo.tbl <- LoadDatabase(inf)
print("Getting genes from database")
start <- Sys.time()
N.sub.lst <- expandingList()
for (jgene in jgenes){
  N.long.filt.query <- filter(motevo.tbl, gene == jgene)  # peaks are not indexed, so dont take them
  N.sub.tmp <- collect(N.long.filt.query, n = Inf)
  N.sub.lst$add(N.sub.tmp)
}
N.long.filt <- N.sub.lst$as.list()
N.long.filt <- bind_rows(N.long.filt)
rm(N.sub.tmp, N.sub.lst)  # worth it? 

# filter peaks after querying database
N.long <- subset(N.long.filt, peak %in% jpeaks & dist <= distfilt)

# Check cooperativityu ----------------------------------------------------

# jgenes <- c("Arntl", "Clock", "Npas2", "Nfil3")

jmotifs <- c("HIC1.p2", "RORA.p2")

N.sub <- subset(N.long, gene %in% jgenes & motif %in% jmotifs) %>%
  group_by(gene, motif, peak) %>%
  summarise(sitecount = sum(sitecount))

# N.sub <- subset(N.long, gene %in% jgenes & motif %in% jmotifs) %>%
#   group_by(gene, motif) %>%
#   summarise(sitecount = sum(sitecount))

# fill missing motifs with 0s
N.sub <- FillWithZeros(N.sub)

N.pairs <- N.sub %>%
  group_by(gene, peak) %>%
  summarise(sitecount = prod(sitecount)) %>%
  arrange(desc(sitecount)) %>%
  group_by(gene) %>%
  summarise(sitecount = sum(sitecount))


# Label amp and phase -----------------------------------------------------

amps.hash <- hash(as.character(fits.long$gene), fits.long$amp.avg)
phases.hash <- hash(as.character(fits.long$gene), fits.long$phase.avg)

N.pairs$amp <- sapply(N.pairs$gene, function(g) amps.hash[[g]])
N.pairs$phase <- sapply(N.pairs$gene, function(g) phases.hash[[g]])


library(PhaseHSV)
jcols <- hsv(PhaseToHsv(2 * pi * N.pairs$phase / 24, 0, 2 *pi), s=0.9, v=0.7)
N.pairs$color <- jcols

pdf(paste0("/home/yeung/projects/tissue-specificity/plots/hic1_rora_cooperativity/hic1_rora_sitecount_vs_amp_phase.disfilt.", distfilt, ".cutoff.", jcutoff, ".pdf"))
ggplot(N.pairs, aes(x = phase, y = amp, size = sitecount)) + geom_point() 
ggplot(N.pairs, aes(x = phase, y = sitecount, size = amp)) + geom_point() 
ggplot(N.pairs, aes(x = sitecount, y = amp, colour = color, label = gene)) + geom_point(size = 2) + geom_text(nudge_y = 0.05) + 
  scale_color_identity() + theme(legend.position = "none")
dev.off()
