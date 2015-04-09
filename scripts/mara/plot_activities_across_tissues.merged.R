# plot_activities_across_tissues.merged.R
# Jake Yeung
# 2015-04-08
# After running cbind_activities.R, we get a matrix of activities across all 
# tissues. Plot polar coordinate figure of the activities. Possibly even do SVD analysis

# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")

IsRnaseq <- function(label.samp){
  # Split by period, if splits into two elements, it's RNASeq, otherwise it's array
  split.length <- length(strsplit(label.samp, "[.]")[[1]])
  if (split.length == 2){
    return(TRUE)
  } else if (split.length == 1){
    return(FALSE)
  }
}

GetMergedColnames <- function(cnames.merged){
  # When we load merged table, our colnames have .1 to represent the RNA-Seq.
  # Fix so that it is Sample.Experiment e.g. Adr18.array or Adr18.rnaseq
  cnames.fixed <- sapply(cnames.merged, function(s){
    # Add .array suffix or .rnaseq suffix depending on if it is Rnaseq or Array
    s.label <- strsplit(s, "[.]")[[1]][[1]]
    if (IsRnaseq(s)){
      s.fixed <- paste0(s.label, ".rnaseq")
    } else {
      s.fixed <- paste0(s.label, ".array")
    }
    return(s.fixed)
  })
}

FitRhythmicWeighted <- function(df, T = 24){
  # Input: long df with exprs, time, se, experiment (array or rnaseq) for
  # a given tissue and a given gene. 
  # 
  # Fits a weighted regression model. With weights in SE
  w = 2 * pi / T  # omega
  tissue <- unique(df$tissue)
  
  sigma.sq <- df$se ^ 2
  jweights <- 1 / sigma.sq
  fit.rhyth <- lm(exprs ~ 0 + experiment + sin(w * time) + cos(w * time), data = df, weights = jweights)
  fit.flat <- lm(exprs ~ 0 + experiment, data = df, weights = jweights)
  ftest <- anova(fit.flat, fit.rhyth)
  ftest.pval <- ftest[["Pr(>F)"]][[2]]
  amp <- unname(sqrt(coef(fit.rhyth)["sin(w * time)"] ^ 2 + coef(fit.rhyth)["cos(w * time)"] ^ 2))
  phase <- unname(atan2(coef(fit.rhyth)["sin(w * time)"], coef(fit.rhyth)["cos(w * time)"] ^ 2))
  df.out <- data.frame(tissue = tissue, amp = amp, phase = phase, pval = ftest.pval)
}

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
ggplot(subset(act.long, gene == jgene), 
       aes(x = time, y = exprs, group = experiment, colour = experiment)) +
  geom_line() +
  geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
  facet_wrap(~tissue) + 
  xlab("CT") +
  ylab("Activity") + 
  ggtitle(jgene)


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

