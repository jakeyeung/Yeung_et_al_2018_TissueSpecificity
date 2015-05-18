# plot_activities_across_tissues.merged.R
# Jake Yeung
# 2015-04-08
# After running cbind_activities.R, we get a matrix of activities across all 
# tissues. Plot polar coordinate figure of the activities. Possibly even do SVD analysis

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
library(dplyr)  # problems with plyr old package
library(reshape2)

# Load data ---------------------------------------------------------------


# merged.act.path <- "results/MARA/MARA_N_centered_with_SE_with_merged/merged/rhythmic_pval1e5amp05/activities.all"
# merged.se.path <- "results/MARA/MARA_N_centered_with_SE_with_merged/merged/rhythmic_pval1e5amp05/standarderrors.all"
# merged.act.path <- "results/MARA/MARA_N_centered_with_SE_with_merged/merged/expressed_genes_threshold5/activities.all"
# merged.se.path <- "results/MARA/MARA_N_centered_with_SE_with_merged/merged/expressed_genes_threshold5/standarderrors.all"
# merged.act.path <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se/rhythmic_pval1e5amp05/activities.all"
# merged.se.path <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se/rhythmic_pval1e5amp05/standarderrors.all"
merged.act.path <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se/expressed_genes_threshold5/activities.all"
merged.se.path <- "/home/yeung/projects/tissue-specificity/results/MARA/MARA_motevo_with_se/expressed_genes_threshold5/standarderrors.all"
merged.act <- read.table(merged.act.path)
merged.se <- read.table(merged.se.path)


# Rename colnames ---------------------------------------------------------

colnames(merged.act) <- GetMergedColnames(colnames(merged.act))
colnames(merged.se) <- GetMergedColnames(colnames(merged.se))

# Create long -------------------------------------------------------------

tissues <- GetTissues.merged(colnames(merged.act))
times <- GetTimes.merged(colnames(merged.act))
experiments <- GetExperiments.merged(colnames(merged.act))

act.long <- data.frame(gene = rep(rownames(merged.act), ncol(merged.act)),
                       tissue = rep(tissues, each = nrow(merged.act)),
                       time = as.numeric(rep(times, each = nrow(merged.act))),
                       exprs = as.numeric(unlist(merged.act)),
                       se = as.numeric(unlist(merged.se)),
                       experiment = rep(experiments, each = nrow(merged.act)))

# save(act.long, file = "Robjs/act.long.Robj")

# Plot  -------------------------------------------------------------------

jgene <- "REST.p3"
jgene <- "HNF1A.p2"
jgene <- "RORA.p2"
jgene <- "SRF.p3"
jgene <- "ARNT_ARNT2_BHLHB2_MAX_MYC_USF1.p2"
jgene <- "AHR_ARNT_ARNT2.p2"
jgene <- "GZF1.p2"
jgene <- "FOX.I1.J2..p2"
jgene <- "ADNP_IRX_SIX_ZHX.p2"
jgene <- "NFIL3.p2"
jgene <- "MEF2.A.B.C.D..p2"
jgene <- "MYOD1.p2"
jgene <- "TBP.p2"

PlotActivitiesWithSE(subset(act.long, gene == jgene))


# Which ones are most rhythmic? -------------------------------------------
act.fit <- act.long %>%
  group_by(tissue, gene) %>%
  do(.data = ., FitRhythmicWeighted(df = .))


# False discovery rate adj ------------------------------------------------

act.fit$pval.adj <- p.adjust(act.fit$pval, method = "BH")


# Show top genes for each tissue ------------------------------------------

head(arrange(data.frame(act.fit), pval), n = 50)


# Plot circle plot across all tissues -------------------------------------

PlotAmpPhaseAllTissues(dat = subset(act.fit, pval.adj <= 5e-2))
PlotAmpPhase(dat = subset(act.fit, pval.adj <= 5e-2 & tissue == "Liver"))

# Plot mean expressions ---------------------------------------------------

act.mean <- act.long %>%
  group_by(tissue, gene, experiment) %>%
  summarise(.data = ., exprs = mean(exprs), se = mean(se))
head(act.mean)

PlotMeanActivitiesWithSE(subset(act.mean, gene == "REST.p3"))


# Do SVD ------------------------------------------------------------------

act.mat <- data.frame(act.long) %>%
  subset(., experiment == "rnaseq") %>%
  dcast(., gene ~ tissue, fun.aggregate = mean, value.var = "exprs")
rownames(act.mat) <- act.mat[, "gene"]
act.mat <- act.mat[, 2:ncol(act.mat)]

act.mat.scale <- t(scale(t(act.mat)))
s <- svd(act.mat)
rownames(s$v) <- colnames(act.mat)
rownames(s$u) <- rownames(act.mat)

plot(s$d^2, type = 'o')

comp <- 1
eigengene <- s$u[, comp]
eigengene2 <- s$u[, comp + 1]
plot(eigengene)
text(eigengene, labels = rownames(names(eigengene)))



# Find most "tissue-specific" ---------------------------------------------

# # use only rnaseq
# act.mean.rnaseq <- subset(act.mean, experiment == "rnaseq")
# tissue.ranges <- ddply(act.mean.rnaseq, .(gene), summarise, range = diff(range(exprs)))
# head(tissue.ranges[order(tissue.ranges$range, decreasing = TRUE), ], n = 10)
# 
# # find liver vs kidney
# act.mean.rnaseq.livkid <- subset(act.mean.rnaseq, tissue %in% c("Liver", "Kidney"))
# tissue.ranges.livkid <- ddply(act.mean.rnaseq.livkid, .(gene), summarise, range = diff(range(exprs)))
# head(tissue.ranges.livkid[order(tissue.ranges.livkid$range, decreasing = TRUE), ], n = 10)
# 
# 
# df <- subset(act.mean, gene == "REST.p3")
# df <- subset(act.mean, gene == "RORA.p2")
# fit.tiss <- lm(exprs ~ 0 + experiment + tissue, data = df)
# fit.null <- lm(exprs ~ 0 + experiment, data = df)
# anova(fit.null, fit.tiss)


# Find rhythmic regulators ------------------------------------------------

# Filter for rhythmic

pval.adj.cutoff <- 0.0005

rhythmic_genes <- act.fit %>%
  group_by(gene) %>%
  mutate(pval.adj.min = min(pval.adj)) %>%
  filter(pval.adj.min < pval.adj.cutoff)

rhythmic_genes <- unique(rhythmic_genes$gene)

act.filt <- subset(act.long, gene %in% rhythmic_genes)

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

act.complex.mat <- dcast(data = act.complex, formula = gene ~ tissue, value.var = "exprs.transformed")
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
for (comp in seq(12)){
  eigengene <- act.svd$v[, comp]
  eigensamp <- act.svd$u[, comp]
  # rotate to phase of largest magnitude in sample of eigengene
  phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
  rotate.factor <- complex(modulus = 1, argument = phase.reference)
  # rotate eigengene by -phase ref
  eigengene <- eigengene * Conj(rotate.factor)
  # rotate eigensamp by +phase ref
  eigensamp <- eigensamp * rotate.factor
  
  v.plot <- PlotComplex2(eigengene, labels = colnames(act.complex.mat), omega = omega, title = paste("Right singular value", comp))
  u.plot <- PlotComplex2(eigensamp, labels = rownames(act.complex.mat), omega = omega, title = paste("Left singular value", comp)) 
  print(v.plot)
  print(u.plot)
}
