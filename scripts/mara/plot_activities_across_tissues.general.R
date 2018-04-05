# plot_activities_across_tissues.general.R
# Jake Yeung
# 2015-06-16
# After running cbind_activities.R, we get a matrix of activities across all 
# tissues. Plot polar coordinate figure of the activities. Possibly even do SVD analysis

# Functions ---------------------------------------------------------------

source("scripts/functions/ActivitiesFunctions.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/SvdFunctions.R")
source("scripts/functions/PlotActivitiesFunctions.R")
library(dplyr)  # problems with plyr old package
library(reshape2)
library(ggplot2)

MaraOutputsToLong <- function(act.path, se.path){
  source("scripts/functions/GetTissueTimes.R")
  act.all <- read.table(act.path)
  se.all <- read.table(se.path)
  tissues <- GetTissues(colnames(act.all))
  times <- GetTimes(colnames(act.all))  
  act.long <- data.frame(gene = rep(rownames(act.all), ncol(act.all)),  # i know it's a motif, bare with me.
                         tissue = rep(tissues, each = nrow(act.all) * length(times)),
                         time = as.numeric(rep(times, length(tissues), each = nrow(act.all))),
                         exprs = as.numeric(unlist(act.all)),
                         se = as.numeric(unlist(se.all)))
}

PlotActivitiesWithSE <- function(act.long, jgene){
  ggplot(data = subset(act.long, gene == jgene) , aes(x = time, y = exprs)) + 
    geom_line() +
    geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
    facet_wrap(~tissue) + 
    xlab("CT") + 
    ylab("Activity") + 
    ggtitle(jgene)
}

# Load data ---------------------------------------------------------------

act.path <- "/home/yeung/projects/tissue-specificity/results/MARA/outputs_MARA_dhs_sitecounts_single_promoters/activities.all"
se.path <- "/home/yeung/projects/tissue-specificity/results/MARA/outputs_MARA_dhs_sitecounts_single_promoters/standarderrors.all"

act.long <- MaraOutputsToLong(act.path, se.path)

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

jgene <- "IRF1.2.7.p3"
jgene <- "ESR1.p3"
jgene <- "LEF1_TCF7_TCF7L1.2.p2"
jgene <- "ZBTB6.p2"

PlotActivitiesWithSE(act.long, jgene)

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
