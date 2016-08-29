# 2016-08-28
# Jake Yeung
# Summarize differences in phases between mRNA and regulator. Include an alpha value?

rm(list=ls())

library(ggplot2)
library(hash)

jmeth <- "g=1001"

source("scripts/functions/ProteomicsFunctions.R")
source("scripts/functions/GetTFs.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

# Load --------------------------------------------------------------------

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.bugfixed.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.long.multimethod.filtbest.staggeredtimepts.bugfixed.annotated.Robj", v=T)
load("Robjs/liver_kidney_atger_nestle/fits.bytiss.bugfixed.Robj", v=T)
fits.long.filt <- subset(fits.long.filt, method == jmeth)
prot.long <- LoadProteomicsData()

prot <- LoadProteomicsData(as.long = FALSE)
prot <- subset(prot, select = c(Gene.names, mean, amp, relamp, phase, pval, qv))
colnames(prot)[colnames(prot) == "Gene.names"] <- "gene"
prot$exper <- "prot"
fits.bytiss$exper <- "rnaseq"

# Get rhythmic TFs and their motifs ---------------------------------------

jmod <- "Liver_SV129,Liver_BmalKO"
genes.liv <- as.character(subset(fits.long.filt, model == jmod)$gene)

tfs <- GetTFs(split.commas = TRUE, get.mat.only = FALSE)

genes.tfs <- genes.liv[genes.liv %in% tfs]
# genes.tfs <- genes.liv

# Can we get TFs that do not have motifs? ---------------------------------
# see if Jingkui have already annotated proteins as TF or not. 


fits.sub <- subset(fits.bytiss, gene %in% genes.tfs & tissue == "Liver_SV129", select = c(gene, amp, phase, pval, exper))
prot.sub <- subset(prot, gene %in% genes.tfs, select = c(gene, amp, phase, pval, exper))
# if gene assigned to two peptides, take one with best rhythmicity
prot.sub <- prot.sub %>%
  group_by(gene) %>%
  filter(pval == min(pval))

merged.long <- rbind(as.data.frame(fits.sub), as.data.frame(prot.sub))
merged <- dcast(merged.long, formula = gene ~ exper, value.var = "phase")

# add worst pvalue
merged.worstpval <- merged.long %>%
  group_by(gene) %>%
  filter(pval == max(pval))
worstpvals <- hash(as.character(merged.worstpval$gene), merged.worstpval$pval)

# add to merged
merged$pval.worst <- sapply(as.character(merged$gene), function(g) worstpvals[[g]])

merged.sub <- subset(merged, pval.worst < 0.05)

ggplot(merged.sub, aes(x = rnaseq, y = prot, label = gene)) +
  geom_text() + theme_bw(18) + geom_abline(slope = 1) + theme(legend.position = "bottom")

ggplot(merged, aes(x = rnaseq, y = prot, label = gene, alpha = -log10(pval.worst))) +
  geom_text() + theme_bw(18) + geom_abline(slope = 1) + theme(legend.position = "bottom")
  

# Add Activity ------------------------------------------------------------


