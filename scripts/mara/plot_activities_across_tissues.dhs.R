# plot_activities_across_tissues.merged.R
# Jake Yeung
# 2015-06-15
# After running cbind_activities.R, we get a matrix of activities across all 
# tissues. Plot polar coordinate figure of the activities. Possibly even do SVD analysis

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/ActivitiesFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
library(dplyr)  # problems with plyr old package
library(reshape2)

# Load data ---------------------------------------------------------------

act.path <- "/home/yeung/projects/tissue-specificity/results/MARA/outputs_MARA_dhs_sitecounts/activities.all"
se.path <- "/home/yeung/projects/tissue-specificity/results/MARA/outputs_MARA_dhs_sitecounts/standarderrors.all"
act.all <- read.table(act.path)
se.all <- read.table(se.path)


# Rename colnames ---------------------------------------------------------

tissues <- GetTissues(colnames(act.all))
times <- GetTimes(colnames(act.all))

act.long <- data.frame(gene = rep(rownames(act.all), ncol(act.all)),  # i know it's a motif, bare with me.
                       tissue = rep(tissues, each = nrow(act.all) * length(times)),
                       time = as.numeric(rep(times, length(tissues), each = nrow(act.all))),
                       exprs = as.numeric(unlist(act.all)),
                       se = as.numeric(unlist(se.all)))


# Plot  -------------------------------------------------------------------

jgene <- "REST.p3"
jgene <- "HNF1A.p2"
jgene <- "HNF$"
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
jgene <- "HNF4A_NR2F1.2.p2"
jgene <- "bHLH_family.p2"
jgene <- "FOXA2.p3"
jgene <- "ATF2.p2"
jgene <- "ZNF423.p2"
jgene <- "HSF1.2.p2"
jgene <- "TGIF1.p2"
jgene <- "ZNF143.p2"
jgene <- "CTCF.p2"
jgene <- "GATA6.p2"

ggplot(data = subset(act.long, gene == jgene) , aes(x = time, y = exprs)) + 
  geom_line() +
  geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
  facet_wrap(~tissue) + 
  xlab("CT") + 
  ylab("Activity") + 
  ggtitle(jgene)


# Which ones are most rhythmic? -------------------------------------------
act.fit <- act.long %>%
  group_by(tissue, gene) %>%
  do(.data = ., FitRhythmicWeighted.singleintercept(df = .))

# False discovery rate adj ------------------------------------------------

act.fit$pval.adj <- p.adjust(act.fit$pval, method = "BH")


# Show top genes for each tissue ------------------------------------------

head(arrange(data.frame(act.fit), pval), n = 50)


# Plot circle plot across all tissues -------------------------------------

PlotAmpPhaseAllTissues(dat = subset(act.fit, pval.adj <= 1))
PlotAmpPhase(dat = subset(act.fit, pval.adj <= 1 & tissue == "Liver"))

# Plot mean expressions ---------------------------------------------------

act.mean <- act.long %>%
  group_by(tissue, gene) %>%
  summarise(.data = ., exprs = mean(exprs), se = mean(se))
head(act.mean)

PlotMeanActivitiesWithSE.singleintercept(subset(act.mean, gene == "REST.p3"))


# Do SVD ------------------------------------------------------------------

act.mat <- data.frame(act.long) %>%
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
eigensamp <- s$v[, comp]
eigensamp2 <- s$v[, comp + 1]

plot(eigengene, eigengene2)
text(eigengene, eigengene2, labels = names(eigengene))

plot(eigensamp, eigensamp2)
text(eigensamp, eigensamp2, labels = names(eigensamp))

