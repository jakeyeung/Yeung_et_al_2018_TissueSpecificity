# 2015-09-07
# Find promoters that are expressed versus not expressed (fitting mixtures of distributions?)

source("scripts/functions/LoadKallisto.R")


# Load --------------------------------------------------------------------

tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")


# Add up by gene ----------------------------------------------------------

tpm.merged.bygene <- tpm.merged %>%
  group_by(gene_name) %>%
  summarise(tpm.gene = sum(tpm))

plot(density(log2(tpm.merged$tpm[which(tpm.merged$tpm > 0)] * 1000 + 1)))
plot(density(log2(tpm.merged.bygene$tpm.gene[which(tpm.merged.bygene$tpm.gene > 0)] * 1000 + 1)))
