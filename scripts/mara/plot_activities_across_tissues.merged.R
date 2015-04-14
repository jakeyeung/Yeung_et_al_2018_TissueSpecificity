# plot_activities_across_tissues.merged.R
# Jake Yeung
# 2015-04-08
# After running cbind_activities.R, we get a matrix of activities across all 
# tissues. Plot polar coordinate figure of the activities. Possibly even do SVD analysis

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesMergedFunctions.R")
source("scripts/functions/PlotFunctions.R")
# Load data ---------------------------------------------------------------


merged.act.path <- "results/MARA/MARA_N_centered_with_SE_with_merged/merged/rhythmic_pval1e5amp05/activities.all"
merged.se.path <- "results/MARA/MARA_N_centered_with_SE_with_merged/merged/rhythmic_pval1e5amp05/standarderrors.all"
merged.act.path <- "results/MARA/MARA_N_centered_with_SE_with_merged/merged/expressed_genes_threshold5/activities.all"
merged.se.path <- "results/MARA/MARA_N_centered_with_SE_with_merged/merged/expressed_genes_threshold5/standarderrors.all"

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

PlotActivitiesWithSE(subset(act.long, gene == jgene))


# Which ones are most rhythmic? -------------------------------------------

act.split <- split(act.long, act.long$tissue)

act.split.fit <- lapply(act.split, function(df.tiss){
  ddply(df.tiss, .(gene), FitRhythmicWeighted)
})


# Combine data ------------------------------------------------------------

act.fit <- do.call(rbind, act.split.fit)

# False discovery rate adj ------------------------------------------------

act.fit$pval.adj <- p.adjust(act.fit$pval, method = "BH")


# Show top genes for each tissue ------------------------------------------

for (t in unique(tissues)){
  print(head(act.split.fit[[t]][order(act.split.fit[[t]]$pval), ]))
}
head(act.fit[order(act.fit$pval.adj), ])


# Plot circle plot across all tissues -------------------------------------

PlotAmpPhaseAllTissues(dat = subset(act.fit, pval.adj <= 5e-2))

# Plot mean expressions ---------------------------------------------------


act.split.mean <- lapply(act.split, function(df.tiss){
  ddply(df.tiss, .(gene), GetExprsMean)
})

act.mean <- do.call(rbind, act.split.mean)

PlotMeanActivitiesWithSE(subset(act.mean, gene == "REST.p3"))


# Plot bar plots across tissues -------------------------------------------

ggplot(theTable,
       aes(x=reorder(Position,Position,
                     function(x)-length(x)))) +
  geom_bar()

ggplot(df.test, 
       aes(x = reorder(gene, exprs), y = exprs)) +
  geom_bar(stat = "identity")

# Find most "tissue-specific" ---------------------------------------------

# use only rnaseq
act.mean.rnaseq <- subset(act.mean, experiment == "rnaseq")
tissue.ranges <- ddply(act.mean.rnaseq, .(gene), summarise, range = diff(range(exprs)))
head(tissue.ranges[order(tissue.ranges$range, decreasing = TRUE), ], n = 10)

# find liver vs kidney
act.mean.rnaseq.livkid <- subset(act.mean.rnaseq, tissue %in% c("Liver", "Kidney"))
tissue.ranges.livkid <- ddply(act.mean.rnaseq.livkid, .(gene), summarise, range = diff(range(exprs)))
head(tissue.ranges.livkid[order(tissue.ranges.livkid$range, decreasing = TRUE), ], n = 10)


df <- subset(act.mean, gene == "REST.p3")
df <- subset(act.mean, gene == "RORA.p2")
fit.tiss <- lm(exprs ~ 0 + experiment + tissue, data = df)
fit.null <- lm(exprs ~ 0 + experiment, data = df)
anova(fit.null, fit.tiss)
