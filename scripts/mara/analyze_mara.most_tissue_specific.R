# Jake Yeung
# analyze_mara.R
# 2015-06-23
setwd("/home/yeung/projects/tissue-specificity")

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

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se.redo/activities"
act.long <- LoadActivitiesLong(indir)


# Which ones are most rhythmic? -------------------------------------------
act.fit <- act.long %>%
  group_by(tissue, gene) %>%
  do(.data = ., FitRhythmicWeighted(dat = ., T = 24, intercepts = TRUE))


# Which are most tissue-specific? -----------------------------------------

act.flat <- act.fit %>%
  group_by(gene) %>%
  summarise(Range = diff(range(rnaseq.int))) %>%
  arrange(desc(Range))


# pdf("plots/activities/tissue_specific_activities.pdf")
# for(m in act.flat$gene){
#   print(PlotActivitiesWithSE(subset(act.long, gene == m)))
# }
# dev.off()


# PCA on regulators -------------------------------------------------------

M <- dcast(act.fit, formula = gene ~ tissue, value.var = "rnaseq.int")
rownames(M) <- M$gene
M$gene <- NULL
M <- t(scale(t(M), center = TRUE, scale = FALSE))
s <- prcomp(M)
screeplot(s)
biplot(s)

# And by rhythms? ---------------------------------------------------------

act.rhyth <- act.fit %>%
  group_by(gene) %>%
  summarise(min.pval = min(pval)) %>%
  arrange(min.pval)

pdf("plots/activities/rhythmic_activities.pdf")
mapply(function(jgene, jtitle) print(PlotActivitiesWithSE(subset(act.long, gene == jgene), jtitle = jtitle)), act.rhyth$gene, paste(act.rhyth$gene, signif(act.rhyth$min.pval, digits=3)))
dev.off()
