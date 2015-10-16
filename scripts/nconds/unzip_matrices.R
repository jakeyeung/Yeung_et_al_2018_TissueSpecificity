setwd("/home/yeung/projects/tissue-specificity")

source("scripts/functions/NcondsFunctions.R")

chunkmain <- "/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_max_3_chunks_20000"
outdir <- "/home/yeung/projects/tissue-specificity/data/nconds2/mats_11_tiss_max_3_chunks_20000.unzipped"
chunks <- list.files(chunkmain)

chunks.unzipped <- UnzipMatrices(file.path(chunkmain, chunk), outdir = outdir)
lapply(chunks, function(chunk){
  UnzipMatrices(file.path(chunkmain, chunk), outdir = outdir)
})

# source("scripts/functions/NcondsFunctions.R")
# mat.chunk.mat <- lapply(mat.chunk, function(des.mat){
#   des.mat$mat <- as.matrix(des.mat$mat)
#   return(des.mat)
# })
start <- Sys.time()
fits <- LoadDesMatDatGeneRunFits(dat.gene, mat.chunk.mat, sparse = FALSE, normalize.weights = FALSE)
print(Sys.time() - start)

start <- Sys.time()
fits <- LoadDesMatDatGeneRunFits(dat.gene, mat.chunk, sparse = TRUE, normalize.weights = FALSE)
print(Sys.time() - start)