LoadActivitiesLong <- function(indir, act.file="activities.all", se.file="standarderrors.all", shorten.motif.name=FALSE){
  source("~/projects/tissue-specificity/scripts/functions/ActivitiesMergedFunctions.R")
  source("~/projects/tissue-specificity/scripts/functions/GetTissueTimes.R")
  source("~/projects/tissue-specificity/scripts/functions/RemoveP2Name.R")
  merged.act <- read.table(file.path(indir, act.file))
  merged.se <- read.table(file.path(indir, se.file))
  
  # Rename colnames ---------------------------------------------------------
  
  colnames(merged.act) <- GetMergedColnames(colnames(merged.act))
  colnames(merged.se) <- GetMergedColnames(colnames(merged.se))
  
  # Create long -------------------------------------------------------------
  
  tissues <- GetTissues.merged(colnames(merged.act))
  times <- GetTimes.merged(colnames(merged.act))
  experiments <- GetExperiments.merged(colnames(merged.act))
  
  if (shorten.motif.name){
    rownames(merged.act) <- sapply(rownames(merged.act), RemoveP2Name)
  }
  
  act.long <- data.frame(gene = rep(rownames(merged.act), ncol(merged.act)),
                         tissue = rep(tissues, each = nrow(merged.act)),
                         time = as.numeric(rep(times, each = nrow(merged.act))),
                         exprs = as.numeric(unlist(merged.act)),
                         se = as.numeric(unlist(merged.se)),
                         experiment = rep(experiments, each = nrow(merged.act)))
  return(act.long)
}

GetTimesTissuesGenoKL <- function(cnames){
  # Get TImes Tissues and Genotypes for Kidney and Liver
  times <- lapply(cnames, function(cname){
    # either BmalKO or SV129
    if (grepl("BmalKO", cname)){
      geno <- "BmalKO"
    } else if (grepl("SV129", cname)){
      geno <- "SV129"
    } else {
      warning("Neither BmalKO or SV129")
    }
    time <- as.numeric(strsplit(cname, geno)[[1]][[2]])
    tissue <- strsplit(cname, paste0("_", geno))[[1]][[1]]
    return(list(time = time, geno = geno, tissue = tissue))
  })
}

LoadActivitiesLongKidneyLiver <- function(indir, act.file="activities.all", se.file="standarderrors.all", collapse.geno.tissue=TRUE, shorten.motif.name=TRUE){
  # handle for Kidney and Liver SV129 and BmalKO
  source("~/projects/tissue-specificity/scripts/functions/ActivitiesMergedFunctions.R")
  source("~/projects/tissue-specificity/scripts/functions/GetTissueTimes.R")
  source("~/projects/tissue-specificity/scripts/functions/RemoveP2Name.R")
  source('scripts/functions/LiverKidneyFunctions.R')
  merged.act <- read.table(file.path(indir, act.file))
  merged.se <- read.table(file.path(indir, se.file))

  if (shorten.motif.name){
    rownames(merged.act) <- sapply(rownames(merged.act), RemoveP2Name)
  }
  time.geno.tissue <- GetTimesTissuesGenoKL(colnames(merged.act))
  tissues <- sapply(time.geno.tissue, function(ll) ll[["tissue"]])
  genos <- sapply(time.geno.tissue, function(ll) ll[["geno"]])
  times <- sapply(time.geno.tissue, function(ll) ll[["time"]])
  
  act.long <- data.frame(gene = rep(rownames(merged.act), ncol(merged.act)),
                         geno = rep(genos, each = nrow(merged.act)),
                         tissue = rep(tissues, each = nrow(merged.act)),
                         time = as.numeric(rep(times, each = nrow(merged.act))),
                         exprs = as.numeric(unlist(merged.act)),
                         experiment = "rnaseq",
                         se = as.numeric(unlist(merged.se)))
  if (collapse.geno.tissue){
    act.long <- CollapseTissueGeno(act.long)
  }
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