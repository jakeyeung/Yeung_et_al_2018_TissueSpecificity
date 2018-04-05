# 2015-11-09
# gaussian_centers.R
# Get Gaussian centers to find alternative promoter usage.

library(ellipse)
# Function ----------------------------------------------------------------

source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotUCSC.R")


# Load --------------------------------------------------------------------
load("Robjs/tpm.afe.avg.Robj", verbose=T)
load("Robjs/dat.long.fixed_rik_genes.Robj")
load("Robjs/tpm.merged.Robj", verbose=T)

load("Robjs/tpm.gauss.Robj", verbose=T)
load("Robjs/tpm.mr.fuzzy.Robj", verbose=T)

# Sort --------------------------------------------------------------------

tpm.gauss.scores <- subset(tpm.gauss, amp.range > 0.5, select = -sigs)
tpm.gauss.scores <- tpm.gauss.scores[order(tpm.gauss.scores$center.dists, decreasing = TRUE), ]
head(tpm.gauss.scores, n = 50)

tpm.mr <- tpm.mr[order(tpm.mr$inter.score, decreasing=TRUE), ]
data.frame(head(tpm.mr, n = 50))


# Look at hits ------------------------------------------------------------

hits <- head(tpm.gauss.scores$gene_name, n = 100)


for (jgene in head(tpm.gauss.scores$gene_name, n = 20)){
  PromoterSpacePlots(tpm.afe.avg, jgene)
}
for (jgene in head(tpm.mr$gene_name, n = 20)){
  PromoterSpacePlots(tpm.afe.avg, jgene)
}
# jgene <- "Ddc"
# PlotGeneAcrossTissues(subset(dat.long, gene == jgene))


# Scale and mean center ---------------------------------------------------

jgene <- "Slc45a3"
PromoterSpacePlots(subset(tpm.afe.avg, mean >= 0), jgene)

PlotGeneAcrossTissues(subset(dat.long, gene == jgene))

dat.sub <- subset(dat.long, gene == "Pvalb")
dat.sub.scaled <- dat.sub %>%
  group_by(tissue) %>%
  mutate(exprs = scale(exprs, center = TRUE, scale = TRUE))
PlotGeneAcrossTissues(dat.sub.scaled)


proms.full <- GetPromoterUsage(subset(tpm.afe.avg, gene_name == jgene), jvar = "tpm_norm.avg", do.svd = T, append.tiss = TRUE, get.means = TRUE, get.entropy = FALSE)  
