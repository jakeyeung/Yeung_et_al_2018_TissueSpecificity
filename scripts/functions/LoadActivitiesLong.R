LoadActivitiesLong <- function(indir, act.file="activities.all", se.file="standarderrors.all"){
  source("~/projects/tissue-specificity/scripts/functions/ActivitiesMergedFunctions.R")
  source("~/projects/tissue-specificity/scripts/functions/GetTissueTimes.R")
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

LoadActivitiesLongDhs <- function(indir, act.file, se.file){
  # expect columns to be just tissues (no time).
  source("~/projects/tissue-specificity/scripts/functions/ActivitiesMergedFunctions.R")
  source("~/projects/tissue-specificity/scripts/functions/GetTissueTimes.R")
  merged.act <- read.table(file.path(indir, act.file))
  merged.se <- read.table(file.path(indir, se.file))
  
  # Create long -------------------------------------------------------------
  
  tissues <- colnames(merged.act)
  act.long <- data.frame(gene = rep(rownames(merged.act), ncol(merged.act)),
                         tissue = rep(tissues, each = nrow(merged.act)),
                         exprs = as.numeric(unlist(merged.act)),
                         se = as.numeric(unlist(merged.se)))
  return(act.long)
}