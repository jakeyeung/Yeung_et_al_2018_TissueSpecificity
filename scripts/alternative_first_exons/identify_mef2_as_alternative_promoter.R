# 2015-07-11
# Find Mef2c as alternative promoter usage?

library(ggplot2)
source("scripts/functions/LoadKallisto.R")
source("scripts/functions/AlternativeFirstExonsFunctions.R")

# Load data ---------------------------------------------------------------

tpm.merged <- LoadKallisto(path.kallisto = "data/alternative_exon_usage/abundance.merged.annotated.sorted.pythonmerged.bed")

# find tissue-specific circadian genes ------------------------------------------------------------

# dat.rhyth <- FitRhythmicDatLong(dat.long)
# dat.rhyth$is.rhythmic <- mapply(IsRhythmic2, pval = dat.rhyth$pval, avg.exprs = dat.rhyth$int.rnaseq)
# dat.rhyth <- dat.rhyth %>%
#   group_by(gene) %>%
#   do(IsTissueSpecific2(., pval.min = 1e-5, pval.max = 0.05, cutoff = 6))
# 
# tissue.spec.circ.genes <- unique(dat.rhyth[which(dat.rhyth$is.tissue.spec.circ == TRUE), ]$gene)

# Calculate isoform usage -------------------------------------------------

tpm.afe <- tpm.merged %>%
  group_by(gene_name, tissue, time) %>%
  mutate(tpm_normalized = CalculateFractionIsoformUsage(tpm, pseudocount = 1))


# Explore mef2 ------------------------------------------------------------

tpm.mef2 <- subset(tpm.afe, gene_name == "Mef2c")

ggplot(tpm.mef2, aes(x = time, y = log2(100 * tpm + 1), group = transcript_id, colour = transcript_id, fill = transcript_id)) + 
  geom_point() + geom_line() + facet_wrap(~tissue)
