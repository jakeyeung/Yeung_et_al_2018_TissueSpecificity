# 2016-07-01
# Jake Yeung
# Alternative promoter usage
# 2016-07-26: redo after merging fastqs (bugfixed)

rm(list=ls())

setwd("/home/yeung/projects/tissue-specificity")

library(dplyr)
library(ggplot2)

source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")

eps <- 1  # for log2 transform

# Load --------------------------------------------------------------------

inf="/home/shared/atgerWTKO_kidneyWTKO/LiverKidney_SV129BmalKO_RF_TotalRNAPolyARNA_kallisto_abundances.annotated.bugfixed.txt"
dat <- read.table(inf, header = TRUE)

genes <- dat$gene_name
transcripts <- dat$target_id
sampnames <- colnames(dat)[grepl("^Kidney|^Liver", colnames(dat))]
dat.exprs <- dat[, sampnames]


# Make long ---------------------------------------------------------------

tissues <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[1]], USE.NAMES = FALSE)
experiments <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[2]], USE.NAMES = FALSE)
times <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[3]], USE.NAMES = FALSE)
genos <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[4]], USE.NAMES = FALSE)
feedings <- sapply(sampnames, function(s){
  jrep <- strsplit(s, "_")[[1]][[5]]
  # remove .A
  jrep <- strsplit(jrep, "\\.")[[1]][[1]]
  return(jrep)
}, USE.NAMES = FALSE)
replicates <- sapply(sampnames, function(s){
  r <- tryCatch({
    rep.letter <- strsplit(s, "\\.")[[1]][[2]]
    rep.numb <- as.numeric(chartr("ABCDE", "12345", rep.letter))
    return(rep.numb)
  }, error = function(e) {
    return(NA)
  })
  return(r)
}, USE.NAMES = FALSE)

# Integrate replicates into ZT time
times.new <- mapply(function(time, jrep){
  time.numb <- as.numeric(strsplit(time, "ZT")[[1]][[2]])
  if (!is.na(jrep)){
    time.numb.withrep <- time.numb + 24 * (jrep - 1)
  } else {
    time.numb.withrep <- time.numb - 48  # if NA, it is Liver, which starts at ZT50, make it start at ZT02 instead
  }
  add.leading.zeros <- FALSE
  if (add.leading.zeros){
    # add leading zeros (OPTIONAL)
    time.numb.withrep <- sprintf("%02d", time.numb.withrep)
    time.new <- paste0("ZT", time.numb.withrep)
  } else {
    time.new <- time.numb.withrep
  }
  return(time.new)
}, times, replicates, USE.NAMES = FALSE)

dat.bytranscript <- data.frame(gene = rep(genes, length(sampnames)),
                               transcript = rep(transcripts, length(sampnames)),
                               tpm = unlist(dat.exprs),
                               tissue = rep(tissues, each = length(genes)),
                               time = rep(times.new, each = length(genes)),
                               geno = rep(genos, each = length(genes)),
                               feeding = rep(feedings, each = length(genes)),
                               experiment = "RNASeq")

# save to file 
save(dat.bytranscript, file = "Robjs/liver_kidney_atger_nestle/dat.bytranscript.bugfixed.annotated.Robj")
