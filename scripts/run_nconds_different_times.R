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

conds <- c('Liver', 'Liver')  # order matters

# filter by conditions
conds1.grep <- paste0(conds[1], collapse = "|")
conds2.grep <- paste0(conds[2], collapse = "|")
dat.filt.array <- dat[, grepl(conds1.grep, colnames(dat))]
dat.filt.seq <- rna.seq.exprs[, grepl(conds2.grep, colnames(rna.seq.exprs))]
colnames(dat.filt.seq) <- sapply(colnames(dat.filt.seq), function(col){
  return(paste0("RNASeq", col))
})
dat.filt <- cbind(dat.filt.array, dat.filt.seq)

start <- Sys.time()
A <- nconds(dat.filt, out.prefix = "plots/nconds_different_times", ncores = 20)
print(Sys.time() - start)

# fit <- nconds(dat.filt, conds, t, out.prefix = "nconds_8_tissues")
# 
# data.6 <- PrepareData(dat, conds, t)
# 
# load("~/projects/select-rhythmic-models/fit_7_conditions.Robj") # fit
# # load("~/projects/select-rhythmic-models/matrices_7_conditions.Robj")  # my_mat
# # load("~/projects/select-rhythmic-models/matrices_2_conditions.Robj")  # my_mat
# 
# dat.with.fit <- InsertFitToMat(fit, data.6)
# 
# col_i.relamps <- grepl("relamp", colnames(dat.with.fit))
# col_i.amps <- grepl("^amp", colnames(dat.with.fit))
# relamps <- apply(dat.with.fit, 1, GetColsApply, col_i=col_i.relamps)
# amps <- apply(dat.with.fit, 1, GetColsApply, col_i = col_i.amps)
# 
# 
# # filter out low relamps
# filtered.genes <- names(relamps[which(amps >= 0.5)])  # amps is tunable
# dat.with.fit.filtered <- dat.with.fit[filtered.genes, ]
# 
# # filter by BICW
# dat.with.fit.filtered <- dat.with.fit[which(dat.with.fit$BICW >= 0.2), ]
# 
# plot_models(dat.with.fit.filtered, 
#             file_path_name = "plots/nconds/7_conds_filtered_02_bicw/7_conds_filtered", 
#             t, 
#             co = length(conds), 
#             conds, 
#             period = 24,
#             test = TRUE)
# 
# # write list of gene names to file for gene enrichment analysis later
# sink("plots/nconds/7_conds_filtered_02_bicw/filtered_genes.txt")
# for (gene in filtered.genes){
#   cat(gene)
#   cat("\n")
# }
# sink()

dev_mode(F)

# # Run on 2 conditions -----------------------------------------------------
# 
# load('~/projects/select-rhythmic-models/data/tissue_array_dat.Robj')
# conds <- c("Liver", "Kidney")
# co <- length(conds)
# t <- rep(seq(18, 64, 2), co)
# grep_conds <- paste0(conds, collapse = "|")
# dat <- dat[, grepl(grep_conds, colnames(dat))]
# dat.df <- data.frame(name = rownames(dat), dat)
# nconds(dat)
# nconds(dat, conds, t, prepare.data = TRUE)
# 
