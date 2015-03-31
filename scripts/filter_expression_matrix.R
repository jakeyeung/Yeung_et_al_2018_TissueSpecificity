# Filter expression matrix to contain only genes that are above a threshold.
# Use RNA-Seq only to determine the genes.

genelst.outpath <- "results/expressed_genes.threshold5.txt"

# Functions ---------------------------------------------------------------

source("scripts/functions/LoadAndHandleData.R")

rowMax <- function(rnaseq){
  apply(rnaseq, 1, max)
}

# Main --------------------------------------------------------------------

rnaseq <- LoadRnaSeq(rna.seq.path = "data/rna_seq_deseq_counts_colnames_fixed.txt")

rnaseq <- log2(rnaseq + 1)

plot(density(unlist(rnaseq)))  # threshold ~5
abline(v = 5)

threshold <- 5  # log2 scale

rnaseq.filt <- rnaseq[which(rowMax(rnaseq) > threshold), ]

# write rownames to file
genes.exprsed <- data.frame(Gene.ID = rownames(rnaseq.filt))

write.table(genes.exprsed, file = genelst.outpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
