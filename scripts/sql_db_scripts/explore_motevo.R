# Explore MOTEVO database
# Jake Yeung
# 

rm(list=ls())
library(dplyr)
source("/home/yeung/projects/sleep_deprivation/scripts/functions/DatabaseFunctions.R")

inf <- "/home/shared/sql_dbs/closestbed_multiple_genes.genomewide.merged.motifindexed.sqlite3"
motevo.tbl <- LoadDatabase(inf)



# Example: Dbp ------------------------------------------------------------

jgene <- "Egr2"
# (explore.query <- filter(atacseq.tbl, !sample %in% c("ZT03_NSD_1", "ZT03_NSD_2", "ZT03_NSD_3", "ZT03_NSD_4")))
(explore.query <- filter(atacseq.tbl, abs(dist) < 1000))  # take all genes near promoters and then normalize !