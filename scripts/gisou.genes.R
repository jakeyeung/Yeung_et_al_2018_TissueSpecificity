# gisou.genes.R
# Are Gisou's genes circadian?
# From email: SwissPalm access


# Source functions --------------------------------------------------------

source(file.path("scripts", "functions", "LoadAndHandleData.R"))
source(file.path("scripts", "functions", "MergeToLong.R"))
source(file.path("scripts", "functions", "PlotGeneAcrossTissues.R"))

# Define dirs -------------------------------------------------------------

# define dirs
data.dir <- "data"
normalized.array.fname <- "array.adj.0.07.txt"
normalized.array.path <- file.path(data.dir, normalized.array.fname)
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)
gisous.genes.path <- file.path(data.dir, "PATs_APTs_uniprot.txt")

# Load file ---------------------------------------------------------------

array <- LoadNormalizedArray(normalized.array.path, remove.negs = TRUE)
rna.seq <- LoadRnaSeq(rna.seq.path)
gisous.table <- read.table(gisous.genes.path, header = TRUE, sep = "\t")
gisous.genes <- gisous.table$Mouse.Homolog


# Log2 transform of array and rnaseq --------------------------------------

normalized.array <- log2(normalized.array + 1)
rna.seq.exprs <- log2(rna.seq.exprs + 1)


# Create data frame -------------------------------------------------------

dat.array.seq <- MergeToLong(normalized.array, rna.seq.exprs)


# Plot genes of interest --------------------------------------------------

dat.array.seq.sub <- subset(dat.array.seq, gene %in% gisous.genes)

pdf("plots/gisou_genes_over_time.pdf")
for (jgene in gisous.genes){
  if (complete.cases(array[jgene, ])){
    print(PlotGeneAcrossTissues(subset(dat.array.seq.sub, gene == jgene), jtitle = jgene)) 
  }
}
dev.off()
