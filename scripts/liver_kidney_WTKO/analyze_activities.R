# 2016-06-27
# Jake Yeung

rm(list=ls())

load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T)

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
source("scripts/functions/LoadActivitiesLong.R")
source("scripts/functions/RemoveP2Name.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source('scripts/functions/LiverKidneyFunctions.R')

library(dplyr)  # problems with plyr old package
library(reshape2)


# Load --------------------------------------------------------------------

# indir <- "/home/yeung/projects/tissue-specificity/results/MARA/atger_kidney/atger_with_kidney"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA/atger_kidney.Liver_SV129/atger_with_kidney.Liver_SV129"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA/atger_kidney.Liver_SV129,Kidney_SV129/atger_with_kidney.Liver_SV129,Kidney_SV129"
# indir <- "/home/yeung/projects/tissue-specificity/results/MARA/atger_kidney.Kidney_SV129/atger_with_kidney.Kidney_SV129"

# with enhancers 40kb
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/atger_kidney.dist_40000.cutoff_3.model_Liver_SV129.method_g=1001/atger_kidney.Liver_SV129.g=1001"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/atger_kidney.dist_40000.cutoff_3.model_Liver_SV129.method_g=1001.cross_TRUE/atger_kidney.Liver_SV129.g=1001"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/atger_kidney.dist_40000.cutoff_3.model_Liver_SV129.method_g=1001.cross_TRUE/atger_kidney.Liver_SV129.g=1001"

act.long <- LoadActivitiesLongKidneyLiver(indir, collapse.geno.tissue = TRUE, shorten.motif.name = TRUE)

jgene <- "RORA"
jgene <- "CUX2"
PlotActivitiesWithSE.rnaseq(subset(act.long, gene == jgene & experiment == "rnaseq")) + theme_bw(24)

# find most tissue-specific
load("Robjs/liver_kidney_atger_nestle/dat.long.liverkidneyWTKO.Robj", v=T); dat.orig <- dat.long
dat.long <- StaggeredTimepointsLivKid(dat.long)
dat.long <- CollapseTissueGeno(dat.long, keep.tissue.col = TRUE)
act.long$tiss <- sapply(as.character(act.long$tissue), function(tiss) strsplit(tiss, "_")[[1]][[1]])

act.tiss <- subset(act.long) %>%
  group_by(tiss, gene) %>%
  summarise(exprs = mean(exprs))

act.diff <- subset(act.tiss) %>%
  group_by(gene) %>%
  summarise(exprs.diff = diff(exprs))
act.diff <- act.diff[order(act.diff$exprs), ]
act.diff$gene <- factor(as.character(act.diff$gene), levels = unique(as.character(act.diff$gene)))
ggplot(act.diff, aes(y = exprs.diff, x = gene, label = gene)) + geom_text() + theme(axis.text = element_blank())

# also true for mRNA?
dat.tiss <- subset(dat.long) %>%
  group_by(gene, tissue.old) %>%
  summarise(exprs = mean(exprs))

source("scripts/functions/GetTFs.R")
tf.genes <- GetTFs()
dat.diff <- subset(dat.tiss) %>%
  group_by(gene) %>%
  summarise(exprs.diff = -diff(exprs))  # positive means higher in liver
dat.diff <- dat.diff[order(dat.diff$exprs.diff), ]
dat.diff.tf <- subset(dat.diff, gene %in% tf.genes)

omega <- 2 * pi / 24
act.complex <- act.long %>%
  group_by(gene, tissue) %>%
  do(ProjectToFrequency2(., omega, add.tissue=TRUE))

s.act <- SvdOnComplex(act.complex, value.var = "exprs.transformed")

jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
for (comp in seq(1)){
  eigens.act <- GetEigens(s.act, period = 24, comp = comp, adj.mag = TRUE, constant.amp = 4, label.n = 45)
  multiplot(eigens.act$v.plot, eigens.act$u.plot, layout = jlayout)
}

print(head(sort(Mod(eigens.act$eigensamp), decreasing=TRUE), n = 50))
jgene <- "NRF1"
jgene <- "ELK1.4_GABP.A.B1."
jgene <- "PAX2"
jgene <- "RORA"
jgene <- "bHLH_family"
jgene <- "SRF"
jgene <- "MAFB"
jgene <- "MTF1"
jgene <- "KLF4"

PlotActivitiesWithSE.rnaseq(subset(act.long, gene == jgene & experiment == "rnaseq")) + theme_bw(24)
