# Create Y input for MARA
# creaet_y_input_for_elastic_net.R
# February 26 2015


# Functions ---------------------------------------------------------------

source("scripts/functions/RemoveExtension.R")
source("scripts/functions/LoadAndHandleData.R")
source("scripts/functions/FixGeneName.R")

WriteHeaders <- function(){
  # write header for elastic input
  cat(paste0("\t", "cos.part", "\t", "sin.part", "\n"))
}


# Get genome wide genes ---------------------------------------------------
array <- LoadNormalizedArray(normalized.array.path = "data/array.adj.0.07.txt", remove.negs = TRUE)
genes <- rownames(array)
genes <- unlist(sapply(genes, FixGeneName))

# Load fit.list from previous code ----------------------------------------

load("/home/yeung/projects/tissue-specificity/results/fits.Robj")  # fit.list from find.oscillating.genes.R


# Write Y outputs ---------------------------------------------------------

conds <- c('Adr', 'Aorta', 'BFAT', 'Kidney', 'Liver', 'Lung', 'Mus')  # each tissue its own Y input

outdir <- "y_input_elastic_net"
outprefix <- "genome_wide/genome_wide"
for (cond in conds){
  sink(file.path(outdir, paste0(outprefix, ".", cond, ".elasticinput")))
  WriteHeaders()
  for (gene in genes){
    amp <- fit.list[[cond]][[gene]]$amp
    phase <- fit.list[[cond]][[gene]]$phase
    if (is.null(phase)){
      next
    }
    # sin.part <- amp * cos(phase)  # if phase comes from model: y = Asin(w + phi)
    # cos.part <- amp * sin(phase)  # if phase comes from y = Asin(w + phi)
    sin.part <- amp * sin(phase)  # # if phase comes from y = Acos(w - phi)
    cos.part <- amp * cos(phase)  # if phase comes from y = Acos(w - phi)
    cat(paste0(gene, "\t", cos.part, "\t", sin.part, "\n"))
  }
  sink()
} 

