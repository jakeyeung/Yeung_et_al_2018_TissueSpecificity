# 2016-02-19
# get_tissue_spec_peaks.R
# Find tissue specific peaks within liver-specific rhythmic genes

rm(list=ls())

library(hash)
library(dplyr)
library(reshape2)

# Functions ---------------------------------------------------------------

source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/DataHandlingFunctions.R")
source("scripts/functions/LdaFunctions.R")

IsTissSpec <- function(is.upper.vec, compare.vec){
  if (all(is.upper.vec == compare.vec)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# Load --------------------------------------------------------------------

load("Robjs/S.long.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)
load("Robjs/N.long.all_genes.all_signif_motifs.Robj", v=T)

# Get liver genes ---------------------------------------------------------

ampmin <- 0.25
jmodels <-c("Kidney;Liver", "Kidney,Liver")
jgenes <- unique(as.character(subset(fits.best, model %in% jmodels & amp.avg > ampmin)$gene))
bg.genes <- unique(as.character(subset(fits.best, model == "")$gene))

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
cutoffs.tiss.lower <- hash(as.character(S.tissuecutoff$tissue), as.numeric(S.tissuecutoff$cutoff))

jstart <- Sys.time()  # 4 min 
S.sub.collapse <- subset(S.long, gene %in% c(jgenes, bg.genes))

# apply is faster than mapply
S.sub.collapse$is.upper <- apply(S.sub.collapse, 1, function(row) IsSignalUpper(as.character(row[5]), as.numeric(row[6]), cutoffs.tiss.upper))
# S.sub.collapse$is.lower <- apply(S.sub.collapse, 1, function(row) IsSignalLower(as.character(row[5]), as.numeric(row[6]), cutoffs.tiss.lower))

S.sub.collapse.peaks <- subset(S.sub.collapse, is.upper == TRUE)
# S.sub.collapse.peaks.lower <- subset(S.sub.collapse, is.lower == TRUE)

compare.vec <- c(F, F, T, T, F, F)

livkid.peaks <- as.character(subset(S.sub.collapse.peaks, tissue %in% c("Liver", "Kidney") & is.upper == TRUE)$peak)
S.sub.collapse.tissspec <- subset(S.sub.collapse, peak %in% livkid.peaks) %>%
  group_by(peak, gene) %>%
  summarise(is.tissspec = IsTissSpec(is.upper, compare.vec)) %>%
  filter(is.tissspec == TRUE)

(nrow(S.sub.collapse))
(nrow(S.sub.collapse.peaks))



# Make sitecounts ---------------------------------------------------------

head(N.long.filt)

N.sub <- subset(N.long.filt, peak %in% S.sub.collapse.tissspec$peak)

# run LDA
mat.fg <- dcast(subset(N.sub, gene %in% jgenes), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bg <- dcast(subset(N.sub, gene %in% bg.genes), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)

# bind rows create labels
mat.fgbg <- bind_rows(mat.fg, mat.bg)
mat.fgbg[is.na(mat.fgbg)] <- 0
rownames(mat.fgbg) <- paste(mat.fgbg$peak, mat.fgbg$gene, sep = ";"); mat.fgbg$peak <- NULL; mat.fgbg$gene <- NULL
labels <- c(rep(1, nrow(mat.fg)), rep(2, nrow(mat.bg)))

out <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  

BoxplotLdaOut(out)
PlotLdaOut(out)

# do crosses
mat.fgbg.cross <- CrossProduct(mat.fgbg, remove.duplicates = TRUE)
dim(mat.fgbg.cross)
# add single factors
mat.fgbg.cross <- cbind(mat.fgbg, mat.fgbg.cross)
# remove columns with 0 variance 
mat.fgbg.cross[which(colSums(mat.fgbg.cross) == 0)] <- list(NULL)

out.cross <- PenalizedLDA(mat.fgbg.cross, labels, lambda = 0.01, K = 1, standardized = FALSE)

BoxplotLdaOut(out.cross)
# PlotLdaOut(out.cross)
m <- SortLda(out.cross)

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
