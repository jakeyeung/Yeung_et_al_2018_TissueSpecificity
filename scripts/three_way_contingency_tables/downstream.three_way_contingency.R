# 2016-10-13
# Jake Yeung

rm(list=ls())

library(dplyr)
library(ggplot2)

# Load --------------------------------------------------------------------

K <- 300
prefix <- "Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.K."
prefix <- "Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K."
inf <- paste0(prefix, K, ".Robj")
inf <- paste0("Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.", K, ".withNmatall.Robj")
load(inf, v=T)

print(fits %>% arrange(pval))

# Which pairs with ONECUT?? -----------------------------------------------

jmotif <- "bHLH_family"
jmotif <- "FOXA"
jmotif <- "RORA"
jmotif <- "CUX2"
jmotif <- "ONECUT"

fits.sub <- subset(fits, grepl(jmotif, pair)) %>% arrange(pval)
print(data.frame(head(fits.sub, n = 15)))

# top guys
head(fits %>% arrange(pval))

# Plot outputs ------------------------------------------------------------

# top guys for FOXA2, ONECUT CUX2
jmotifs <- c("FOXA2", "ONECUT1.2", "CUX2")
grepstr <- paste(jmotifs, collapse = "|")

(fits.sub <- subset(fits, grepl(grepstr, pair)) %>% arrange(pval))

# contains RORA?
fits.sub$has.ROR <- sapply(fits.sub$pair, function(p) grepl("RORA", p))

fits.sub$pair <- factor(fits.sub$pair, levels = fits.sub$pair)
m <- ggplot(subset(fits.sub, pval < 0.05), aes(x = pair, y = -log10(pval), fill = has.ROR)) + geom_bar(stat = "identity") +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Top motif pairs between\nFOXA2, ONECUT, or CUX2") + 
  xlab("")
m 


# Check to see if there are bad pairs -------------------------------------

Npath <- "/home/yeung/projects/tissue-specificity/data/sitecounts/motevo/sitecount_matrix_geneids"
N <- read.table(Npath, header=TRUE)
motifs <- colnames(N)[which(colnames(N) != "Gene.ID")]
source("/home/yeung/projects/tissue-specificity/scripts/functions/GetTFs.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/HandleMotifNames.R")
motifs <- sapply(as.character(motifs), function(m) RemoveP2AndCommaBracesDashes(m))
motifs.output <- colnames(N.mat.all)[!colnames(N.mat.all) %in% c("gene", "peak", "model")]
motifs.common <- intersect(motifs, motifs.output)

motif.pairs <- combn(as.character(motifs.common), m = 2, simplify = FALSE)

N.mat.freqs <- lapply(motif.pairs, function(motif.pair){
  motif.pair.str <- paste(motif.pair, collapse = ";")
  N.mat.freq <- N.mat.all %>%
    group_by_(.dots = c("model", motif.pair)) %>%
    summarise(freq = length(gene)) %>%
    mutate(pair = motif.pair.str)
  return(N.mat.freq)
  # should be a matrix of 8 rows, otherwise something is wrong and you should skip it
  # if (nrow(N.mat.freq) != 8){
  #   # print(paste("Warning: some zeros in sitecounts encountered: skipping", motif.pair.str))
  #   warning(paste("Warning: some zeros in sitecounts encountered: skipping", motif.pair.str))
  #   return(data.frame())
  # } else {
  #   return(N.mat.freq)
  # }
})

library(data.table)
N.mat.freqs.rbind <- rbindlist(N.mat.freqs)
colnames(N.mat.freqs.rbind) <- c("model", "motif1", "motif2", "freq", "pair")

pairs.sub <- as.character(fits.sub$pair)

N.mat.freqs.sum <- subset(N.mat.freqs.rbind, pair %in% pairs.sub) %>%
  group_by(pair) %>%
  summarise(nrows=length(pair)) %>%
  arrange(nrows)

bad.pairs <- subset(N.mat.freqs.sum, nrows != 8)$pair


# Redo plot but remove bad pairs ------------------------------------------

m.fixed <- ggplot(subset(fits.sub, pval < 0.05 & !pair %in% bad.pairs), aes(x = pair, y = -log10(pval), fill = has.ROR)) + geom_bar(stat = "identity") +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Top motif pairs between\nFOXA2, ONECUT, or CUX2") + 
  xlab("")
m.fixed



# Document this change ----------------------------------------------------

save(fits, N.mat.freqs.rbind, N.mat.all, file = "Robjs/three_way_cooccurence/three.way.cooccurrence.bugfixed.nmodels.2.K.300.withNmatallandNmatfreqs.Robj")


# Further exploration -----------------------------------------------------



# how many pairs with FOXA, ONECUT, or CUX2??
fits.sub$partner <- sapply(as.character(fits.sub$pair), function(p){
  pvec <- strsplit(p, ";")[[1]]
  partner <- pvec[!pvec %in% jmotifs]
  if (length(partner) == 0) partner <- NA
  return(partner)
})

fits.sum <- subset(fits.sub, pval < 0.05) %>%
  group_by(partner) %>%
  summarise(n.hits = length(partner)) %>%
  arrange(desc(n.hits))

fits.sum$partner <- factor(fits.sum$partner, levels = fits.sum$partner)
m.sum <- ggplot(subset(fits.sum), aes(x = partner, y = n.hits)) + geom_bar(stat = "identity")

