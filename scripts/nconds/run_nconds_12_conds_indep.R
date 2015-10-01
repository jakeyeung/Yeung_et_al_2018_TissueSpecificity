# Run nconds.R
library("devtools")
dev_mode()
install("~/projects/ncond")  # use jake branch
# load_all("~/projects/ncond")
library(nconds)

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


# Load and plot on 12 conditions -------------------------------------------

# load from this because it is easier, it is run from a code that takes a while to calculate
load("~/projects/select-rhythmic-models/data/tissue_array_dat.Robj", verbose = T)
# conds <- c('Liver', 'BFAT', 'Adr', 'Mus', 'Lung', 'Aorta', 'Kidney')
conds <- c('Adr', 'Aorta', 'BFAT', 'BS', 'Cere', 'Heart', 'Hypo', 'Kidney', 'Liver', 'Lung', 'Mus', 'WFAT')  # order matters
# conds <- c('Adr', 'Aorta', 'BFAT', 'Heart', 'Kidney')
co <- length(conds)
t <- rep(seq(18, 64, 2), co)

# filter by conditions
conds.grep <- paste0(conds, collapse = "|")
dat.filt <- dat[, grepl(conds.grep, colnames(dat))]

start <- Sys.time()
A <- nconds(dat.filt, out.prefix = "~/projects/tissue-specificity/nconds_outputs/nconds_12_tissues/nconds_12_tissues", ncores = 20, write.intermediates = TRUE, all.models = FALSE)
print(Sys.time() - start)

# # make matrix
# period <- 24
# conds <- rep(conds, each = length(unique(t)))
# n.co <- co
# mat <- creat_matrix_list(t, conds, n.co, period, all.models = FALSE)



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
