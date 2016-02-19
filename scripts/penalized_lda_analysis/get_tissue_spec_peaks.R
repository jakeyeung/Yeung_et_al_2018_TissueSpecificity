# 2016-02-19
# get_tissue_spec_peaks.R
# Find tissue specific peaks within liver-specific rhythmic genes

library(hash)
library(dplyr)

# Functions ---------------------------------------------------------------

source("scripts/functions/SitecountsFunctions.R")

# Load --------------------------------------------------------------------

load("Robjs/S.long.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)
load("Robjs/N.long.liver_genes.all_motifs.100000.Robj", verbose=T)  # huge file 

# Get liver genes ---------------------------------------------------------

ampmin <- 0.25
liv.genes <- unique(as.character(subset(fits.best, model == "Liver" & amp.avg > ampmin)$gene))

# shift by 1 log2 unit the fit looks better
log.shift <- 2.5
S.tissuecutoff$cutoff.adj <- 2^(log2(S.tissuecutoff$cutoff) + log.shift)

cutoff.adj <- S.tissuecutoff$cutoff.adj
cutoff.lower <- S.tissuecutoff$cutoff

# show upper and lower limit
pseudo <- 1e-3
ggplot(S.long[sample(x = 1:nrow(S.long), size = 0.01 * nrow(S.long)), ], aes(x = log2(signal + pseudo))) + geom_density() + 
  facet_wrap(~tissue) +
  geom_vline(aes(xintercept = log2(cutoff.adj)), data = S.tissuecutoff, colour = "blue") + 
  geom_vline(aes(xintercept = log2(cutoff)), data = S.tissuecutoff, colour = "red")

cutoffs.tiss.upper <- hash(as.character(S.tissuecutoff$tissue), as.numeric(S.tissuecutoff$cutoff.adj))

jstart <- Sys.time()  # 4 min 
S.sub.collapse <- subset(S.long, gene %in% liv.genes)

# apply is faster than mapply
S.sub.collapse$is.upper <- apply(S.sub.collapse, 1, function(row) IsSignalUpper(as.character(row[5]), as.numeric(row[6]), cutoffs.tiss.upper))

S.sub.collapse.peaks <- subset(S.sub.collapse, is.upper == TRUE)

(nrow(S.sub.collapse))
(nrow(S.sub.collapse.peaks))


# Filter N.long for peaks  ------------------------------------------------

# need to do it tissue by tissue
tissues <- as.character(unique(S.sub.collapse.peaks$tissue))

peaks.bytissue <- lapply(tissues, function(tiss){
  unique(as.character(subset(S.sub.collapse.peaks, tissue == tiss)$peak))
})
names(peaks.bytissue) <- tissues

N.long.sum.bytiss <- lapply(tissues, function(tiss){
  dat.out <- subset(N.long.liver_dhs, peak %in% peaks.bytissue[[tiss]]) %>%
    group_by(gene, motif) %>%
    summarise(sitecount.sum = sum(sitecount))
  dat.out$tissue <- tiss
  return(dat.out)
})
N.long.sum.bytiss <- do.call(rbind, N.long.sum.bytiss)

# normalize by number of peaks for each gene
S.npeaks <- S.sub.collapse.peaks %>%
  group_by(gene, tissue) %>%
  summarise(n.peaks = length(peak))

n.peaks.hash <- hash(paste(S.npeaks$gene, S.npeaks$tissue, sep = ";"), S.npeaks$n.peaks)

# normalize
N.long.sum.bytiss$sitecount.norm.bypeak <- apply(N.long.sum.bytiss, 1, function(row){
  jgene <- as.character(row[1]); tiss <- as.character(row[4]); sitecount <- as.numeric(row[3])
  n.peaks <- n.peaks.hash[[paste(jgene, tiss, sep = ";")]]
  if (is.null(n.peaks)){
    return(NA)
  } else {
    sitecount / n.peaks
  }
})

# remove NAs
N.long.sum.bytiss <- N.long.sum.bytiss[which(!is.na(N.long.sum.bytiss$sitecount.norm.bypeak)), ]

