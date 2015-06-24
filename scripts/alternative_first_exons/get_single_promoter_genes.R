# Jake Yeung
# 2015-06-15
# get_single_promoter_genes.R

library(dplyr)

# functions ---------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/LoadKallisto.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")

IsSinglePromoter <- function(dat){
  # Expect row: tpm_normalized
  # if row all equal, then it suggests identical promoter usage across all tissues
  if (length(unique(dat$transcript_id)) == 1){
    return(data.frame(is.single = TRUE))
  } else {
    return(data.frame(is.single = FALSE))
  }  
}

# Load files --------------------------------------------------------------

load(file = "Robjs/arrayrnaseq.mydat.Robj")
tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")

# Fit Rhythmic ------------------------------------------------------------

dat.rhyth <- FitRhythmicDatLong(mydat)


# Calculate fractional isoform usage --------------------------------------

tpm.afe <- tpm.merged %>%
  group_by(gene_name, tissue, time) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 0)) %>%
  filter(!is.nan(tpm_normalized))  # NaNs come from sum(tpm) == 0, we don't zero-expressed genes

# Which genes only have a single promoter? --------------------------------

# slow ~ 5 minutes
tpm.single.promoters <- tpm.afe %>%
  group_by(gene_name) %>%
  do(IsSinglePromoter(.)) %>%
  filter(is.single == TRUE)

save(tpm.single.promoters, file = "Robjs/tpm.single.promoters.Robj")

tpm.afe.singles <- subset(tpm.afe, gene_name %in% tpm.single.promoters$gene_name)


# Print out my gene list --------------------------------------------------

sink("data/gene_lists/genes_single_promoters")
for (gene in tpm.afe.singles$gene_name){
  cat(gene)
  cat("\n")
}
sink()
  
