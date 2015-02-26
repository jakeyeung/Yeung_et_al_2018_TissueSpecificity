# Create Y input for MARA
# creaet_y_input_for_elastic_net.R
# February 26 2015


# Functions ---------------------------------------------------------------

WriteHeaders <- function(){
  # write header for elastic input
  cat(paste0("\t", "cos.part", "\t", "sin.part", "\n"))
}

# Get genes from Cedric's model -------------------------------------------

genes.fname <- "plots/nconds/7_conds_filtered_02_bicw/7_conds_filtered26.txt"
genes <- ReadListToVector(genes.fname)


# Load fit.list from previous code ----------------------------------------

load("/home/yeung/projects/tissue-specificity/results/fits.Robj")  # fit.list from find.oscillating.genes.R


# Write Y outputs ---------------------------------------------------------

conds <- c('Adr', 'Aorta', 'BFAT', 'Kidney', 'Liver', 'Lung', 'Mus')  # each tissue its own Y input

outdir <- "y_input_elastic_net"
outprefix <- RemoveExtension(basename(genes.fname))
for (cond in conds){
  sink(file.path(outdir, paste0(outprefix, cond, ".elasticinput")))
  WriteHeaders()
  for (gene in genes){
    amp <- fit.list[[cond]][[gene]]$amp
    phase <- fit.list[[cond]][[gene]]$phase
    cos.part <- amp * cos(phase)
    sin.part <- amp * sin(phase)
    cat(paste0(gene, "\t", cos.part, "\t", sin.part, "\n"))
  }
  sink()
}


