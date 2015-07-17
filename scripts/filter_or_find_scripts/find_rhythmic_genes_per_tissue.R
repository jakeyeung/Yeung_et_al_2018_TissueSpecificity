# find_rhythmic_genes_per_tissue.R
# 2015-05-05
# Run linear model to find genes that are rhythmic in one tissue


# Functions ---------------------------------------------------------------

library(hash)
library(dplyr)
source("scripts/functions/LoadArrayRnaSeq.R")
source("scripts/functions/FitRhythmic.R")
source("scripts/functions/GetClockGenes.R")

WriteGeneListToFile <- function(genes, outfile, tissue, ampphasedic){
  sink(file = outfile)
  for (g in genes){
    if (missing(ampphasedic) & missing(tissue)){
      cat(g)
      cat("\n")
    } else{
      jkey = paste(tissue, g, sep = ';')
      jval <- ampphasedic[[jkey]]
      amp <- strsplit(jval, ';')[[1]][[1]]
      phase <- strsplit(jval, ';')[[1]][[2]]
      pval <- strsplit(jval, ';')[[1]][[3]]
      cat(g, "\t", amp, "\t", phase, "\t", pval, "\t")
      cat('\n')
    }
  }
  sink()
}

# Load data ---------------------------------------------------------------


dat.long <- LoadArrayRnaSeq()

# dat.test <- subset(dat.long, gene %in% GetClockGenes())


# Find rhythmic genes -----------------------------------------------------


start <- Sys.time()
dat.fit <- dat.long %>%
  group_by(tissue, gene) %>%
  split(.$tissue) %>%
  mclapply(., function(jdf) do(.data = jdf, FitRhythmic(.)), mc.cores = 12) %>%
  do.call(rbind, .)
print(Sys.time() - start)


# Create hash to relate tissue;gene with amp and phase --------------------

jkeys <- paste(dat.fit$tissue, dat.fit$gene, sep = ";")
jvals <- paste(dat.fit$amp, dat.fit$phase, dat.fit$pval, sep = ';')
ampphasedic <- hash(keys = jkeys, values = jvals)

# Filter NaNs --------------------------------------------------------------


# filter
dat.fit.filt <- subset(dat.fit, !is.na(pval))

plot(density(dat.fit$pval[which(dat.fit$pval > 0)]))
ggplot(data = subset(dat.fit, pval < 1e-5 & amp > 0.5), aes(x = phase)) + geom_density()


# Find top hits -----------------------------------------------------------


# Count how many pass a threshold
pval.cutoff <- 1e-4

top.hits <- dat.fit %>%
  group_by(tissue) %>%
  mutate(pval.adj = p.adjust(pval)) %>%
  split(.$tissue) %>%
  lapply(., function(dat) dat[order(dat$pval), ])

pval.cutoff <- 0.05
n.hits <- lapply(top.hits, function(dat){
  length(which(dat$pval.adj < pval.cutoff))
})
barplot(as.vector(unlist(n.hits)), names.arg = names(n.hits), main = pval.cutoff)
abline(h = 200)


# Take top hits OR top 200 genes ------------------------------------------

take.top.N <- 200
for (tiss in names(top.hits)){
#   outname <- paste0(tiss, "_rhythmic_genes_", pval.cutoff)
  outname <- paste0(tiss, "_rhythmic_genes_with_ampphase", pval.cutoff)
  outfile <- file.path("results/rhythmic_genes_by_tissue", outname)
  if (n.hits[[tiss]] < take.top.N){
    # Return top 200 genes
    rhythmic.genes <- top.hits[[tiss]]$gene[1:take.top.N]
  } else {
    rhythmic.genes <- top.hits[[tiss]]$gene[which(top.hits[[tiss]]$pval.adj < pval.cutoff)]
  }
  print(paste0("Writing ", length(rhythmic.genes), " genes to ", outfile))
  WriteGeneListToFile(rhythmic.genes, outfile, tissue = tiss, ampphasedic = ampphasedic)
}


