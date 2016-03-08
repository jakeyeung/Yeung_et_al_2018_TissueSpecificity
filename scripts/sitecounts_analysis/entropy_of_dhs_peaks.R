# 2016-03-04
# Are RORA elements in Liver-rhythmic genes showing signs of DHS peaks?

library(ggplot2)
library(hash)

source("scripts/functions/ShannonEntropy.R")


# Function ----------------------------------------------------------------


LiverToOthersRatio <- function(jsignal, indx = 4){
  return(jsignal[4] / mean(jsignal[-4]))
}

# Load --------------------------------------------------------------------

load("Robjs/N.long.all_genes.all_signif_motifs.Robj", v=T)
load("Robjs/S.long.Robj", verbose=T)
load("Robjs/fits.best.max_3.collapsed_models.amp_cutoff_0.15.phase_sd_maxdiff_avg.Robj", verbose=T)


# Get liver-rhythmic genes ------------------------------------------------

minamp <- 0.5
jmodel <- c("Kidney;Liver", "Kidney,Liver")
minamp <- 0.25

jmodel <- c("Adr;Liver")

jmodel <- "Liver"
genes.list <- as.character(subset(fits.best, model == jmodel & amp.avg > minamp)$gene)


# Sort peaks by RORA ------------------------------------------------------

jmotif <- "ONECUT1,2.p2"
jmotif <- "MYFfamily.p2"
jmotif <- "ZNF148.p2"
jmotif <- "FOXA2.p3"
jmotif <- "ESR1.p3"
jmotif <- "FOXA2.p3"
jmotif <- "RXRG_dimer.p3"
jmotif <- "HNF4A_NR2F1,2.p2" 
jmotif <- "RORA.p2"
jmotif <- "bHLH_family.p2"

N.sub <- subset(N.long.filt, gene %in% genes.list & motif == jmotif)
N.sub <- N.sub[order(N.sub$sitecount, decreasing=TRUE), ]


# Correlation between entropy and ROR sitecounts --------------------------

jpeaks <- unique(N.sub$peak)

S.entropy <- subset(S.long, peak %in% jpeaks) %>%
  group_by(peak, gene) %>%
  summarise(H=ShannonEntropy(p = signal, normalize = TRUE),
            LRatio = LiverToOthersRatio(jsignal = signal, indx = 4),
            L = signal[4])

entropy.hash <- hash(as.character(S.entropy$peak), as.numeric(S.entropy$H))
LR.hash <- hash(as.character(S.entropy$peak), as.numeric(S.entropy$LRatio))
L.hash <- hash(as.character(S.entropy$peak), as.numeric(S.entropy$L))

N.sub$H <- sapply(as.character(N.sub$peak), function(p) entropy.hash[[p]])
N.sub$LR <- sapply(as.character(N.sub$peak), function(p) LR.hash[[p]])
N.sub$L <- sapply(as.character(N.sub$peak), function(p) L.hash[[p]])

# ggplot(N.sub, aes(x = sitecount, y = H)) + geom_point(alpha = 0.1) + ggtitle(paste("Entropy vs sitecount for", jmotif))
# ggplot(N.sub, aes(x = H)) + geom_histogram() + ggtitle(paste("Entropy histogram for", jmotif))
# ggplot(N.sub, aes(x = LR)) + geom_histogram() + ggtitle(paste("Entropy histogram for", jmotif))

# ggplot(N.sub, aes(x = LR, y = H)) + geom_point(alpha = 0.1) + ggtitle(paste("Entropy vs Liver signal to others ratio", jmotif))
ggplot(N.sub, aes(x = sitecount, y = L)) + geom_point(alpha = 0.1) + ggtitle(paste("Motif vs Liver signal", jmotif))


