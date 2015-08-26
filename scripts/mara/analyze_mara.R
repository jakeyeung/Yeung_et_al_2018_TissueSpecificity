# Jake Yeung
# analyze_mara.R
# 2015-06-23

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")

library(dplyr)  # problems with plyr old package
library(reshape2)

# Load --------------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto_motevo"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/outputs_MARA_dhs_sitecounts_array_kallisto"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1.dist_filt"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto.threshold-5.572458.promoters/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1.union.dist_filt"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto.threshold-5.572458.dist_filt/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto.threshold-5.572458.promoters/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/tissue_specific_rhythmic_genes.pvalmax0.05pvalmin1e-04relamp0.1.dist_filt/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/tissue_specific_rhythmic_genes.pvalmax0.05pvalmin1e-05relamp0.1.dist_filt/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/Liver_rhythmic_genes.pvalmax0.05pvalmin1e-05.relamp0.1"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/Liver_rhythmic_genes.pvalmax0.05pvalmin1e-05relamp0.1"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/Kidney_rhythmic_genes.pvalmax0.05pvalmin1e-04relamp0.1"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/Lung_rhythmic_genes.pvalmax0.05pvalmin0.001relamp0.1"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/Heart_rhythmic_genes.pvalmax0.05pvalmin0.001relamp0.1"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto-threshold-5.572458"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1.union.promoters/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto.threshold-5.572458.promoters"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_deseq"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_deseq"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_deseq_int"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long <- LoadActivitiesLong(indir)


# Plot stuff --------------------------------------------------------------


PlotActivitiesWithSE(subset(act.long, gene == "HNF1A.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "ZNF423.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "ZNF384.p2"))

PlotActivitiesWithSE(subset(act.long, gene == "ATF5_CREB3.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "NFE2L2.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "HAND1.2.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "SOX17.p2"))

PlotActivitiesWithSE(subset(act.long, gene == "HBP1_HMGB_SSRP1_UBTF.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "ZBTB16.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "ALX1.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "FOXA2.p3"))
PlotActivitiesWithSE(subset(act.long, gene == "NR4A2.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "TEAD1.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "RFX1..5_RFXANK_RFXAP.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "MYFfamily.p2"))

PlotActivitiesWithSE(subset(act.long, gene == "RORA.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "LEF1_TCF7_TCF7L1.2.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "HNF1A.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "REST.p3"))
PlotActivitiesWithSE(subset(act.long, gene == "MEF2.A.B.C.D..p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "CTCF.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "HNF4A_NR2F1.2.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "HNF4A_NR2F1.2.p2" & tissue %in% c("Liver", "Kidney")))
PlotActivitiesWithSE(subset(act.long, gene == "ONECUT1.2.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "SOX17.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "NANOG.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "ATF2.p2"))
# PlotActivitiesWithSE(subset(act.long, gene == "ATF6.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "ELK1.4_GABP.A.B1..p3"))
PlotActivitiesWithSE(subset(act.long, gene == "FOXA2.p3"))
PlotActivitiesWithSE(subset(act.long, gene == "MYOD1.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "NFIL3.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "FOXD3.p2"))


# SVD on RNASeq motif activities ------------------------------------------

act.rnaseq <- act.long %>%
  subset(., experiment == "rnaseq") %>% 
  group_by(tissue, gene) %>%
  summarise(exprs.avg = mean(exprs))

act.rnaseq.mat <- dcast(act.rnaseq, formula = gene ~ tissue, value.var = "exprs.avg")
rownames(act.rnaseq.mat) <- act.rnaseq.mat$gene
act.rnaseq.mat$gene <- NULL

act.svd <- prcomp(t(scale(t(act.rnaseq.mat), center = TRUE, scale = FALSE)), center = FALSE, scale. = FALSE)

screeplot(act.svd, type = "lines")
biplot(act.svd, choices = c(1,2))

# Which ones are most rhythmic? -------------------------------------------
act.fit <- act.long %>%
  group_by(tissue, gene) %>%
  do(FitRhythmicWeighted(dat = .))

# False discovery rate adj ------------------------------------------------

act.fit$pval.adj <- p.adjust(act.fit$pval, method = "BH")

# Show top genes for each tissue ------------------------------------------

head(arrange(data.frame(act.fit), pval), n = 50)

# Plot circle plot across all tissues -------------------------------------

pval.adj.cutoff <- 1

PlotAmpPhaseAllTissues(dat = subset(act.fit, pval.adj <= pval.adj.cutoff))
PlotAmpPhase(dat = subset(act.fit, pval.adj <= pval.adj.cutoff& tissue == "Heart"))
PlotAmpPhase(dat = subset(act.fit, pval.adj <= pval.adj.cutoff& tissue %in% c("BFAT", "Mus"))) + facet_wrap(~tissue)
PlotAmpPhase(dat = subset(act.fit, pval.adj <= pval.adj.cutoff& tissue == "Liver"))
PlotAmpPhase(dat = subset(act.fit, pval.adj <= pval.adj.cutoff& tissue == "Aorta"))

# Find rhythmic regulators ------------------------------------------------


rhythmic_genes <- act.fit %>%
  group_by(gene) %>%
  mutate(pval.adj.min = min(pval.adj)) %>%
  filter(pval.adj.min < pval.adj.cutoff)

rhythmic_genes <- unique(rhythmic_genes$gene)

act.filt <- subset(act.long, gene %in% rhythmic_genes)
PlotAmpPhaseAllTissues(dat = subset(act.fit, gene %in% rhythmic_genes))

library(plyr)
source("scripts/functions/SvdFunctions.R")  # many script-specific functions here

omega <- 2 * pi / 24

start.time <- Sys.time()
act.complex <- lapply(split(act.filt, act.filt$tissue), function(x){
  ddply(x, .(gene), ProjectToFrequency2, omega = omega, add.tissue = TRUE)
}) %>%
  do.call(rbind, .) %>%
  mutate(magnitude = Mod(exprs.transformed)) %>%
  arrange(desc(magnitude))
print(Sys.time() - start.time)
# head(act.complex, n = 100)

detach("package:plyr", unload=TRUE)

library(dplyr)


# SVD on complex matrix ---------------------------------------------------
# optinally filter out motifs

filter.motifs <- c("HNF4A_NR2F1.2.p2", "HNF1A.p2")
# act.complex.mat <- dcast(data = act.complex, formula = gene ~ tissue, value.var = "exprs.transformed")
act.complex.mat <- dcast(data = subset(act.complex, !gene %in% filter.motifs), formula = gene ~ tissue, value.var = "exprs.transformed")
rownames(act.complex.mat) <- act.complex.mat[, "gene"]
act.complex.mat <- act.complex.mat[, 2:ncol(act.complex.mat)]

act.svd <- svd(act.complex.mat) 

# add row and colnames
rownames(act.svd$u) <- rownames(act.complex.mat)
rownames(act.svd$v) <- colnames(act.complex.mat)


# Plot first component ----------------------------------------------------


# screeplot
plot(act.svd$d ^ 2 / sum(act.svd$d ^ 2), type = 'o')  # eigenvalues

theme_set(theme_gray())
# for (comp in seq(1, ncol(act.complex.mat))){
for (comp in seq(5)){
  eigengene <- act.svd$v[, comp]
  eigensamp <- act.svd$u[, comp]
  # rotate to phase of largest magnitude in sample of eigengene
  phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
  rotate.factor <- complex(modulus = 1, argument = phase.reference)
  # rotate eigengene by -phase ref
  eigengene <- eigengene * Conj(rotate.factor)
  # rotate eigensamp by +phase ref
  eigensamp <- eigensamp * Conj(rotate.factor)
  
  v.plot <- PlotComplex2(eigengene, labels = colnames(act.complex.mat), omega = omega, title = paste("Right singular value", comp))
  u.plot <- PlotComplex2(eigensamp, labels = rownames(act.complex.mat), omega = omega, title = paste("Left singular value", comp)) 
  print(v.plot)
  print(u.plot)
}

