# Functions for handling Atger-Kidney from Nestle and Jerome

CollapseTissueGeno <- function(dat.long, keep.tissue.col=FALSE){
  # Collapse tissue into tissuegeno, so you can use functions that expect "tissues" as conditions in colname (hack)
  if (keep.tissue.col){
    dat.long$tissue.old <- dat.long$tissue
  }
  dat.long$tissue <- factor(paste(dat.long$tissue, dat.long$geno, sep = "_"), levels = c("Liver_SV129", "Kidney_SV129", "Liver_BmalKO", "Kidney_BmalKO"))
  return(dat.long)
}

SameTimepointsLivKid <- function(dat.long){
  # Make Liver and Kidney same number of timepoints (tissue column). Should work for collapse or not collapse
  conds <- unique(as.character(dat.long$tissue))
  conds.lst = list()
  for (cond in conds){
    conds.lst[[cond]] <- unique(subset(dat.long, tissue == cond)$time)
  }
  common.times <- Reduce(intersect, conds.lst)
  return(subset(dat.long, time %in% common.times))
}

# Functions for handling DAT_I_I liver kidney Cedric
TissueFromCname <- function(cname){
  # "X0_1 -> Kidney"
  # "X0_1_liver -> Liver"
  jsplit <- strsplit(cname, "_")[[1]]
  if (jsplit[[length(jsplit)]] == "liver"){
    tiss <- "Liver"
  } else {
    tiss <- "Kidney"
  }
  return(tiss)
} 

TimeFromCname <- function(cname, rm.first.char=TRUE){
  # X0_1 -> 0 + (1 - 1) * 24
  # X22_2 -> 22 + (2 - 1) * 24
  if (rm.first.char){
    cname <- substr(cname, start = 2, stop = nchar(cname))
  }
  time.base <- as.numeric(strsplit(cname, "_")[[1]][[1]])
  time.multiplier <- as.numeric(strsplit(cname, "_")[[1]][[2]]) - 1
  time <- time.base + time.multiplier * 24
  return(time)
}

LoadLivKid <- function(inf){
  if (missing(inf)){
    inf <- "/home/yeung/projects/tissue-specificity/data/gene_exprs/liver_v_kidney/DAT_l_l.txt"
  }
  dat.mat <- read.table(inf, header = TRUE)
  genes <- dat.mat$name
  exprs <- subset(dat.mat, select = c(-name))
  tissues <- rep(sapply(colnames(exprs),TissueFromCname, USE.NAMES = FALSE), each = nrow(exprs))
  tissues.uniq <- unique(tissues)
  times <- rep(sapply(colnames(exprs), TimeFromCname, USE.NAMES = FALSE), each = nrow(exprs))
  
  # make long
  dat <- data.frame(gene = dat.mat$name, 
                    tissue = tissues,
                    time = times,
                    exprs = unlist(exprs),
                    experiment = "rnaseq")
}