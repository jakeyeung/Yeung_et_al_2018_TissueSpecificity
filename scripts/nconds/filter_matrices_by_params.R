# 2015-10-07
# I ran 11 tissues max 4 params to generate matrices.
# I have alrdy run lm fit for 11 tissues max 3 params.
# Filter out matrices from 11 tiss max 4 to contain only
# matrices containing rhyth param of 4.

n.params <- 4
chunk.size <- 200000
chunk.count <- 1
chunk.i <- 1
setwd("~/projects/tissue-specificity")

start <- Sys.time()

# Source ------------------------------------------------------------------

source("scripts/functions/ListFunctions.R")

# Load --------------------------------------------------------------------

matdir <- "/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_max_4_chunks_2e+05"
outdir <- "/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_only_4_chunks_2e+05.filtered"
dir.create(outdir)

chunks.paths <- list.files(matdir) 

des.mats.list.filt <- expandingList()
for (chunk.path in chunks.paths){
  print(paste("Loading", file.path(matdir, chunk.path)))
  print(Sys.time() - start)
  load(file.path(matdir, chunk.path))  # des.mats.list
  print("loaded")
  print(Sys.time() - start)
  for (des.mat in des.mats.list){
    n.rhyth <- length(des.mat$rhyth.tiss) - 1  # ignore flat model
    if (n.rhyth == n.params){
      des.mats.list.filt$add(des.mat)
      chunk.count <- chunk.count + 1
    }
    if (chunk.count %% chunk.size == 0){
      # write to file
      fname <- paste0("chunk", chunk.i, ".filt.Robj")
      des.mats.list <- des.mats.list.filt$as.list()
      print(paste("Saving to:", file.path(outdir, fname)))
      save(des.mats.list, file = file.path(outdir, fname))
      rm(des.mats.list.filt.flat, des.mats.list.filt)
      des.mats.list.filt <- expandingList()  # renew
      chunk.i <- chunk.i + 1
    }
  }
}