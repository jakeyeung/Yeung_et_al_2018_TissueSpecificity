# 2015-10-05
# Jake Yeung
# load_chunks_run_nconds.R


set.seed(0)
w <- 2 * pi / 24

library(dplyr)
library(hash)
library(Matrix)

setwd("~/projects/tissue-specificity")


# Source ------------------------------------------------------------------

source("scripts/functions/GetClockGenes.R")
source("scripts/functions/NcondsFunctions.R")
source("scripts/functions/Queue.R")
source("scripts/functions/ListFunctions.R")
source("scripts/functions/PlotGeneAcrossTissues.R")
source("scripts/functions/MemoryManagement.R")
source("scripts/functions/EnvironmentsFunctions.R")


# Objs dir ----------------------------------------------------------------

chunks.path <- "/home/yeung/projects/tissue-specificity/data/datlong_chunks_by_gene/expressed_genes"

dat.long.path <- "Robjs/dat.long.fixed_rik_genes.Robj"

