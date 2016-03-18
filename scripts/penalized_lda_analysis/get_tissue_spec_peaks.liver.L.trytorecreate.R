# 2016-02-19
# get_tissue_spec_peaks.R
# Find tissue specific peaks within liver-specific rhythmic genes
# but compare KL vs Others (use rhythmic genes always)

rm(list=ls())

library(hash)
library(dplyr)
library(reshape2)

# Functions ---------------------------------------------------------------

source("scripts/functions/SitecountsFunctions.R")
source("scripts/functions/DataHandlingFunctions.R")
source("scripts/functions/LdaFunctions.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

# Load --------------------------------------------------------------------

load("Robjs/dat.long.fixed_rik_genes.Robj", v=T)
load("Robjs/S.long.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)
load("Robjs/S.tissuecutoff.Robj", verbose=T)
load("Robjs/N.long.all_genes.all_signif_motifs.Robj", v=T)

# Get liver genes ---------------------------------------------------------

ampmin <- 0.25
jmodels <-c("Liver")
compare.vec <- c(F, F, F, T, F, F)  # Cere, Heart, Kidney, Liver, Lung, Mus
compare.vec.all <- c(T, T, T, T, T, T)
jgenes <- unique(as.character(subset(fits.best, model %in% jmodels & amp.avg > ampmin)$gene))
compare.vec.other <- !compare.vec

# shift by 1 log2 unit the fit looks better
log.shift <- 2.5
# log.shift <- 1
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
S.sub.collapse <- subset(S.long, gene %in% c(jgenes))

# apply is faster than mapply
S.sub.collapse$is.upper <- apply(S.sub.collapse, 1, function(row) IsSignalUpper(as.character(row[5]), as.numeric(row[6]), cutoffs.tiss.upper))

# all peaks with high signal
S.sub.collapse.peaks <- subset(S.sub.collapse, is.upper == TRUE)

# get livkid peaks within foreground genes and other peaks from foreground genes
livkid.peaks <- unique(as.character(subset(S.sub.collapse.peaks, gene %in% jgenes & tissue %in% c("Liver") & is.upper == TRUE)$peak))
bg.peaks <- unique(as.character(subset(S.sub.collapse.peaks, gene %in% jgenes & ! tissue %in% c("Liver") & is.upper == TRUE)$peak))

S.sub.collapse.livkid <- subset(S.sub.collapse, peak %in% livkid.peaks) %>%
  group_by(peak, gene) %>%
#   summarise(is.tissspec = IsTissSpec(is.upper, compare.vec, compare.vec2)) %>%
  summarise(is.tissspec = IsTissSpec(is.upper, compare.vec)) %>%
  filter(is.tissspec == TRUE)

S.sub.collapse.bg <- subset(S.sub.collapse, peak %in% bg.peaks) %>%
  group_by(peak, gene) %>%
  summarise(is.tissspec = IsTissSpec(is.upper, compare.vec.all, reverse = FALSE)) %>%  # take peaks present across all tissues
#   summarise(is.tissspec = IsTissSpec(is.upper, compare.vec, reverse = TRUE)) %>%
  filter(is.tissspec == TRUE)

print(paste("N foreground peaks:", length(unique(as.character(S.sub.collapse.livkid$peak)))))
print(paste("N background peaks:", length(unique(as.character(S.sub.collapse.bg$peak)))))

# coming from how many genes?
print(paste("Number of genes considered:", length(unique(S.sub.collapse.livkid$gene))))
print(unique(S.sub.collapse.livkid$gene))


# Sanity check these peaks ------------------------------------------------

jsub <- subset(S.sub.collapse, peak %in% as.character(S.sub.collapse.livkid$peak) & tissue == "Liver")
plot(density(jsub$signal))
hist(jsub$signal, breaks = 100)


# Do further cuts ---------------------------------------------------------

# keep high signal only
fg.peaks.filt <- as.character(subset(jsub, signal > 2)$peak)
# fg.peaks.filt <- as.character(subset(jsub, signal > 0)$peak)
print(length(fg.peaks.filt))
S.sub.collapse.livkid <- subset(S.sub.collapse.livkid, peak %in% fg.peaks.filt)

# down sample bg.peaks
set.seed(0)
S.sub.collapse.bg <- S.sub.collapse.bg[sample(seq(length(fg.peaks.filt)), size = length(fg.peaks.filt)), ]

print("After further cutoffs")
print(paste("N foreground peaks:", length(unique(as.character(S.sub.collapse.livkid$peak)))))
print(paste("N background peaks:", length(unique(as.character(S.sub.collapse.bg$peak)))))

# show one randomluy
sample(unique(as.character(S.sub.collapse.bg$peak)), 1)
sample(unique(as.character(S.sub.collapse.livkid$peak)), 1)


# Sanity check background genes -------------------------------------------

# do they have high signal in kidney?
jsub.bg <- subset(S.sub.collapse, peak %in% as.character(S.sub.collapse.bg$peak) & tissue == "Kidney")

plot(density(jsub.bg$signal))
hist(jsub.bg$signal, breaks = 100)

# Make sitecounts ---------------------------------------------------------

head(N.long.filt)
N.sub.livkid <- subset(N.long.filt, peak %in% S.sub.collapse.livkid$peak)
N.sub.other <- subset(N.long.filt, peak %in% S.sub.collapse.bg$peak)

# run LDA
mat.fg <- dcast(subset(N.sub.livkid), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)
mat.bg <- dcast(subset(N.sub.other), formula = peak + gene ~ motif, value.var = "sitecount", fun.aggregate = sum, fill = 0)

# bind rows create labels
mat.fgbg <- bind_rows(mat.fg, mat.bg)
mat.fgbg[is.na(mat.fgbg)] <- 0
rownames(mat.fgbg) <- paste(mat.fgbg$peak, mat.fgbg$gene, sep = ";"); mat.fgbg$peak <- NULL; mat.fgbg$gene <- NULL
labels <- c(rep(1, nrow(mat.fg)), rep(2, nrow(mat.bg)))

out <- PenalizedLDA(mat.fgbg, labels, lambda = 0.1, K = 1, standardized = FALSE)  

BoxplotLdaOut(out, jtitle = "Single factor separation")
PlotLdaOut(out)
m.singles <- SortLda(out)
# do crosses
colnames(mat.fgbg) <- sapply(colnames(mat.fgbg), RemoveP2Name)
mat.fgbg.cross <- CrossProduct(mat.fgbg, remove.duplicates = TRUE)
dim(mat.fgbg.cross)
# add single factors
mat.fgbg.cross <- cbind(mat.fgbg, mat.fgbg.cross)
# remove columns with 0 variance 
mat.fgbg.cross[which(colSums(mat.fgbg.cross) == 0)] <- list(NULL)

# jlambda <- 0.029  # kidliv
jlambda <- 0.015  # liv only
out.cross <- PenalizedLDA(mat.fgbg.cross, labels, lambda = jlambda, K = 1, standardized = FALSE)
m <- SortLda(out.cross)
print(length(m))
BoxplotLdaOut(out.cross, jtitle = "Cross product separation")
PlotLdaOut(out.cross, take.n = 50, from.bottom = TRUE)


# Write peaks to output  -------------------------------------------------

