# create_x_input_for_elastic_net.R
# February 26 2014


# Functions ---------------------------------------------------------------

setwd("~/projects/tissue-specificity")
source("scripts/functions/RemoveExtension.R")
source("scripts/functions/FixGeneName.R")
source("scripts/functions/ReadListToVector.R")

# Load site counts and genes ----------------------------------------------

site.mat.full <- read.table("site_count_matrices/sitecounts_matrix_EPD", header = TRUE, row.names = 1)

genes.fname <- "plots/nconds/7_conds_filtered_02_bicw/7_conds_filtered4.txt"
genes <- ReadListToVector(genes.fname)
genes <- unlist(sapply(genes, FixGeneName))

# Filter sites to match genes ---------------------------------------------

site.mat <- site.mat.full[which(rownames(site.mat.full) %in% genes), ]

# Write to file
out.fname <- RemoveExtension(basename(genes.fname))
out.dir <- "site_count_matrices/"
write.table(site.mat, 
            file = file.path(out.dir, paste0(out.fname, ".sitecounts")),
            sep = "\t",
            quote = FALSE,
            col.names = NA)
