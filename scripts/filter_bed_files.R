# FilterBedFiles.R
#
# Given list of genes, filter bed file so it only contains gene of interest.

source("scripts/functions/ReadListToVector.R")
source("scripts/functions/FixGeneName.R")

bedfile <- read.table("bedfiles/Mm_EPDnew.liftOver.mm10.bed")  # 4th column contains gene names V4

genes <- ReadListToVector("plots/nconds/7_conds_filtered_05_amp/7_conds_filtered4.txt")

genes.fixed <- unlist(sapply(genes, FixGeneName))

genes.grep <- paste0(genes.fixed, collapse = "|")

bedfile.filtered <- bedfile[grepl(genes.grep, bedfile$V4), ]

write.table(bedfile.filtered, 
            file = "bedfiles/7_conds_filtered4.bed", 
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
