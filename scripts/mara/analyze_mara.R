# Jake Yeung
# analyze_mara.R
# 2015-06-23

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
library(dplyr)  # problems with plyr old package
library(reshape2)

LoadActivitiesLong <- function(indir, act.file="activities.all", se.file="standarderrors.all"){
  merged.act <- read.table(file.path(indir, act.file))
  merged.se <- read.table(file.path(indir, se.file))
  
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
  return(act.long)
}

# Load --------------------------------------------------------------------

indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto_motevo"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/expressed_genes_kallisto/"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/rhythmic_genes_by_tissues_kallistoarray.pval0.001.relamp0.1"
indir <- "/home/yeung/projects/tissue-specificity/results/MARA/outputs_MARA_dhs_sitecounts_array_kallisto"

act.long <- LoadActivitiesLong(indir)


# Plot stuff --------------------------------------------------------------

PlotActivitiesWithSE(subset(act.long, gene == "RORA.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "RORA.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "LEF1_TCF7_TCF7L1.2.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "HNF1A.p2"))
PlotActivitiesWithSE(subset(act.long, gene == "REST.p3"))
PlotActivitiesWithSE(subset(act.long, gene == "MEF2.A.B.C.D..p2"))
PlotActivitiesWithSE(subset(act.long, gene == "CTCF.p2"))


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

# Find rhythmic regulators ------------------------------------------------

# Filter for rhythmic

pval.adj.cutoff <- 0.05

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
for (comp in seq(length(act.svd$d))){
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
