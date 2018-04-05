# Jake Yeung
# read_dhs_signal_save_to_robj.R
# Read DHS signal (matrix from pipeline /home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/motevo_dhs_scripts_clean from VITALIT)
# 2016-04-05

library(dplyr)

start <- Sys.time()

args <- commandArgs(trailingOnly=TRUE)

inf <- args[1]
outf <- args[2]
pseudo <- 0.01

if (file.exists(outf)) stop("outf exists, not overwriting...")

print("Reading table")
S <- read.table(inf)

# Add colnames ------------------------------------------------------------

tissues <- c("Cere", "Heart", "Kidney", "Liver", "Lung", "Mus")
cnames <- c("chromo", "start", "end", tissues, "gene", "dist")
colnames(S) <- cnames

# Normalized dat ----------------------------------------------------------

for (tiss in tissues){
  S[[tiss]] <- 10^6 * S[[tiss]] / sum(S[[tiss]])
}

print("Mat to long")
signal.vec <- unlist(S[, colnames(S) %in% tissues])
S.long <- data.frame(chromo = S$chromo, start = S$start, end = S$end, 
                     peak = paste(paste(S$chromo, S$start, sep = ":"), S$end, sep = "-"), # chr1:7800234-7800734
                     tissue = rep(tissues, each = nrow(S)), 
                     signal = signal.vec,
                     gene = S$gene, dist = S$dist)

print(paste("Add Z-score with pseudo", pseudo))
S.long <- S.long %>%
  group_by(tissue) %>%
  mutate(zscore = scale(log2(signal + pseudo), center = TRUE, scale = TRUE))

print(paste("Saving to", outf))
save(S.long, file=outf)

print(Sys.time() - start)
