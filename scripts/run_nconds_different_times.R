# Run run_nconds_different_times.R
# Run for different time intervals to see what it looks like
library("devtools")
dev_mode()
install("~/projects/ncond")  # use jake branch
library(nconds)

setwd("~/projects/tissue-specificity")
scripts.dir <- "scripts"
funcs.dir <- "functions"
data.dir <- "data"
source(file.path(scripts.dir, funcs.dir, "LoadAndHandleData.R"))
rna.seq.fname <- "rna_seq_deseq_counts_colnames_fixed.txt"
rna.seq.path <- file.path(data.dir, rna.seq.fname)

GetRelamp <- function(dat.with.fit, jgene){
  # Inputs
  # dat.with.fit: output of InsertFitToMat. Expect colnames "relamp_i"
  # jgene: row to extract
  #
  # Outputs:
  # Relative amps across conditions
  
  relamps <- dat.with.fit[jgene, grepl("relamp_", colnames(dat.with.fit))]
  return(max(relamps))
}

GetColsApply <- function(dat.with.fit.slice, col_i){
  # Inputs
  # dat.with.fit: output of InsertFitToMat. Expect colnames "relamp_i"
  #
  # Outputs:
  # Relative amps across condition
  cols <- dat.with.fit.slice[col_i]
  return(max(cols))
}

AddUnderscores <- function(cnames){
  # rnaseq colnames are Adr22. Change them to Adr_22
  cnames.new <- sapply(cnames, function(s){
    time <- substr(s, nchar(s) - 1, nchar(s))
    tissue <- substr(s, 1, nchar(s) - 2)
    return(paste0(tissue, "_", time))
  })
}

# Load and plot on 2 conditions -------------------------------------------

# load from this because it is easier, it is run from a code that takes a while to calculate
load("~/projects/select-rhythmic-models/data/tissue_array_dat.Robj")  # dat
rna.seq.exprs <- LoadRnaSeq(rna.seq.path)

colnames(rna.seq.exprs) <- AddUnderscores(colnames(rna.seq.exprs))

# get common genes
common.genes <- intersect(rownames(dat), rownames(rna.seq.exprs))
dat <- dat[common.genes, ]
rna.seq.exprs <- rna.seq.exprs[common.genes, ]

conds <- c('Liver')  # order matters

# filter by conditions
conds1.grep <- paste0(conds[1], collapse = "|")
dat.filt.array <- dat[, grepl(conds1.grep, colnames(dat))]

# subset of array as second time point
t.vec <- seq(22, 64, 6)
subset.grep <- paste0(t.vec, collapse = "|")
cnames.subset <- paste0("reviL_", t.vec)

# dat.filt.array2 <- dat.filt.array[, grepl(subset.grep, colnames(dat.filt.array))]
colnames(dat.filt.array2) <- cnames.subset
dat.filt <- cbind(dat.filt.array, dat.filt.array2)

start <- Sys.time()
A <- nconds(dat.filt, out.prefix = "plots/nconds_different_times/nconds_different_times2", ncores = 20)
print(Sys.time() - start)


dev_mode(F)


