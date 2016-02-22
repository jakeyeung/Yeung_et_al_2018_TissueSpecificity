# 2016-02-19
# get_tissue_spec_peaks.R
# Find tissue specific peaks within liver-specific rhythmic genes

library(hash)
library(dplyr)

# Functions ---------------------------------------------------------------

source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/DataHandlingFunctions.R")


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
    summarise(sitecount.sum = sum(sitecount), sitecount.max = max(sitecount))
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
gene.i <- which(colnames(N.long.sum.bytiss) == "gene")
tiss.i <- which(colnames(N.long.sum.bytiss) == "tissue")
sitecount.i <- which(colnames(N.long.sum.bytiss) == "sitecount.sum")

N.long.sum.bytiss$sitecount.norm.bypeak <- apply(N.long.sum.bytiss, 1, function(row){
  jgene <- as.character(row[gene.i]); tiss <- as.character(row[tiss.i]); sitecount <- as.numeric(row[sitecount.i])
  n.peaks <- n.peaks.hash[[paste(jgene, tiss, sep = ";")]]
  if (is.null(n.peaks)){
    return(NA)
  } else {
    sitecount / n.peaks
  }
})

# remove NAs
N.long.sum.bytiss <- N.long.sum.bytiss[which(!is.na(N.long.sum.bytiss$sitecount.norm.bypeak)), ]


# Check distributions of motifs -------------------------------------------

library(ggplot2)

jmotif <- "HNF1A.p2"
jmotif <- "ONECUT1,2.p2"
jmotif <- "CUX2.p2"
jmotif <- "RXRG_dimer.p3"

# The palette with grey:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(subset(N.long.sum.bytiss, motif == jmotif), aes(x = sitecount.norm.bypeak, fill = tissue)) + geom_density(alpha = 0.5) + ggtitle(jmotif) + scale_fill_manual(values=cbPalette)
ggplot(subset(N.long.sum.bytiss, motif == jmotif), aes(x = sitecount.max, fill = tissue)) + geom_density(alpha = 0.5) + ggtitle(jmotif) + scale_fill_manual(values=cbPalette)
ggplot(subset(N.long.sum.bytiss, motif == jmotif & tissue == "Liver"), aes(x = sitecount.max, fill = tissue)) + geom_density(alpha = 0.5) + ggtitle(jmotif) + scale_fill_manual(values=cbPalette)


# Do cross product --------------------------------------------------------

# jmotif <- "RXRG_dimer.p3"
# 
# sitecount.motif <- hash(as.character(subset(N.long.liver_dhs, motif == jmotif)$peak), subset(N.long.liver_dhs, motif == jmotif)$sitecount)
# 
# peak.i <- GetRowIndx(N.long.liver_dhs, "peak")
# sc.i <- GetRowIndx(N.long.liver_dhs, "sitecount")
# 
# # 3 minutes?
# N.long.liver_dhs$sitecount.cross <- apply(N.long.liver_dhs, 1, function(row){
#   sc <- as.numeric(row[sc.i])
#   peak <- row[peak.i]
#   motif.sc <- sitecount.motif[[peak]]
#   if (is.null(motif.sc)) motif.sc <- 0
#   return(motif.sc * sc)
# })
# 
# N.long.sum.bytiss <- lapply(tissues, function(tiss){
#   dat.out <- subset(N.long.liver_dhs, peak %in% peaks.bytissue[[tiss]]) %>%
#     group_by(gene, motif) %>%
#     summarise(sitecount.sum = sum(sitecount), sitecount.max = max(sitecount), sitecount.cross.max = max(sitecount.cross))
#   dat.out$tissue <- tiss
#   return(dat.out)
# })
# N.long.sum.bytiss <- do.call(rbind, N.long.sum.bytiss)
