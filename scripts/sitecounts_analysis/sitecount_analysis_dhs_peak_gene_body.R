# 2015-12-08
# Jake Yeung
# site count analysis on gene body to do motif enrichment

library(ggplot2)

# Load --------------------------------------------------------------------


start <- Sys.time()
# N <- read.table("data/sitecounts/motevo_by_peaks_dhs_gene_bodies/merged.closest.bed", nrows = 10)  # 30 GB
S <- read.table("/home/yeung/data/tissue_specificity/motevo_dhs/dhs_signal/dhs_signal_windows500.chr.sorted.closest.mat")
print(Sys.time() - start)

load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj")

# Add colnames ------------------------------------------------------------

tissues <- c("Cere", "Heart", "Kidney", "Liver", "Lung", "Mus")
cnames <- c("chromo", "start", "end", tissues, "chromo.gene", "start.gene", "end.gene", "gene", "blank", "strand", "dist")
colnames(S) <- cnames


# Normalized dat ----------------------------------------------------------

for (tiss in tissues){
  S[[tiss]] <- 10^6 * S[[tiss]] / sum(S[[tiss]])
}

signal.vec <- unlist(S[, colnames(S) %in% tissues])
S.long <- data.frame(chromo = S$chromo, start = S$start, end = S$end, 
                     peak = paste(paste(S$chromo, S$start, sep = ":"), S$end, sep = "-"), # chr1:7800234-7800734
                     tissue = rep(tissues, each = nrow(S)), 
                     signal = signal.vec, 
                     gene = S$gene, dist = S$dist)


# Regions of interest -----------------------------------------------------

liver.genes <- subset(fits.best, model == "Liver")$gene

S.sub <- subset(S.long, gene %in% liver.genes)

# check for bias across tissues
ggplot(S.long[sample(x = 1:nrow(S.long), size = 0.01 * nrow(S.long)), ], aes(x = log10(signal))) + geom_density() + facet_wrap(~tissue)

# lets look at Celsr1 near chr15:85,959,964-85,964,941
jgene <- "Celsr1"
startmin <- 85959964
endmin <- 85964941

S.subsub <- subset(S.sub, gene == jgene & start > startmin & end < endmin)
