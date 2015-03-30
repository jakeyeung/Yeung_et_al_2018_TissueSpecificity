# plot_activities_across_tissues.R
# Jake Yeung
# 2015-03-30
# After running cbind_activities.R, we get a matrix of activities across all 
# tissues. Plot polar coordinate figure of the activities. Possibly even do SVD analysis


# Functions ---------------------------------------------------------------

source("scripts/functions/GetTissueTimes.R")
source("scripts/functions/SvdFunctions.R")

ProjectToFrequency2 <- function(df, omega){
  # simpler than ProjectToFrequency().
  # expect df to be gene i in tissue c with column names time and exprs
  exprs.transformed <- DoFourier(df$exprs, df$time, omega = omega)
  return(data.frame(exprs.transformed = exprs.transformed))
}

# Load file ---------------------------------------------------------------


# out.fpath from cbind_activities.R
in.fpath <- "/home/yeung/projects/tissue-specificity/results/outputs_all_genes_MARA.swissregulon.corrected/MARA.33motifs.corrected/activities.all"

act.all <- read.table(in.fpath)



# Create long dataframe ---------------------------------------------------

tissues <- GetTissues(colnames(act.all))
times <- GetTimes(colnames(act.all))
act.long <- data.frame(gene = rep(rownames(act.all), ncol(act.all)),  # i know it's a motif, bare with me.
                       tissue = rep(tissues, each = nrow(act.all) * length(times)),
                       time = as.numeric(rep(times, length(tissues), each = nrow(act.all))),
                       exprs = as.numeric(unlist(act.all)))


# Convert to Fourier ------------------------------------------------------

act.split <- split(act.long, act.long$tissue)

omega <- 2 * pi / 24
act.split.proj <- lapply(act.split, function(df){
  ddply(df, .(gene), ProjectToFrequency2, omega = omega)
})
