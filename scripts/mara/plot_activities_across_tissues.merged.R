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

# Load data ---------------------------------------------------------------


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

ggplot(data = subset(act.long, gene == jgene) , aes(x = time, y = exprs)) + 
  geom_line() +
  geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
  geom_errorbar(limits, width=0.25) + 
  facet_wrap(~tissue) + 
  xlab("CT") + 
  ylab("Activity") + 
  ggtitle(jgene)

jgene <- "REST.p3"
jgene <- "HNF1A.p2"
jgene <- "RORA.p2"
ggplot(subset(act.long, gene == jgene), 
       aes(x = time, y = exprs, group = experiment, colour = experiment)) +
  geom_line() +
  geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
  facet_wrap(~tissue) + 
  xlab("CT") +
  ylab("Activity") + 
  ggtitle(jgene)
  
  
