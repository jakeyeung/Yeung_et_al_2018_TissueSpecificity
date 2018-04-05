# Convert33Motifs.R
# Jake Yeung
# 2015-03-20
# get full site count, filter to jonathan's 33 motifs, output it as a new sitecounts

setwd("~/projects/tissue-specificity")
source("scripts/functions/BiomartFunctions.R")
source("scripts/functions/Filter33Motifs.R")

# Main --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
sitecounts.fpath <- args[1]
sitecounts.filtered.out <- args[2]

# sitecounts.fpath <- "/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/sitecount_matrix.full"
# sitecounts.fpath <- "/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/sitecount_matrix.full.genenames"
# sitecounts.filtered.out <- "/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/sitecount_matrix.full.filtered33motifs"
# load("/home/yeung/projects/tissue-specificity/site_count_matrices/GLM/Wang/pm_WT.Robj") # pm_WT.collapsed


sitecounts.full <- read.table(file = sitecounts.fpath, header = TRUE, sep = "\t")
genenames.uniq <- make.names(sitecounts.full$Gene.ID, unique = TRUE)
rownames(sitecounts.full) <- genenames.uniq
sitecounts.full$Gene.ID <- NULL

sitecounts.full.filt <- Filter33Motifs(sitecounts.full)

Gene.ID <- rownames(sitecounts.full.filt)

sitecounts.full.filt <- cbind(Gene.ID, sitecounts.full.filt)

write.table(sitecounts.full.filt, quote = FALSE, sep = "\t", row.names = FALSE, 
            file = sitecounts.filtered.out)
